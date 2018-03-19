#include <Rcpp.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
using namespace Rcpp;
using namespace std;

S4 _C_dtm(std::vector<std::vector<std::string> >& xx, std::vector<int>& ngram,
          bool tf, bool idf) {

    int n = xx.size();
    vector<vector<string> >* xp = &xx;
    if (!(ngram[0] == 1 && ngram[1] == 2)) {
        int from = ngram[0];
        int to   = ngram[1];
        xp = new vector<vector<string> >(n);
        vector<vector<string> >& x = *xp;
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j + from <= xx[i].size(); ++j) {
                stringbuf sb;
                for (int k = j; k < j + from; ++k)
                    sb.sputn(xx[i][k].c_str(), xx[i][k].size());
                x[i].push_back(sb.str()); //todo
                for (int len = from + 1; len < to; ++len) {
                    if (j + len > xx[i].size()) break;
                    sb.sputn(xx[i][j + len - 1].c_str(), xx[i][j + len - 1].size());
                    x[i].push_back(sb.str()); //todo
                }
            }
        }
    }
    vector<vector<string> >& x = *xp;
    vector<vector<int> > widx;
    map<string, int> code;
    vector<string> word;
    IntegerVector dl(n);
    if (!tf) {
        for (int i = 0; i < n; ++i) {
            sort(x[i].begin(), x[i].end());
            auto t = std::unique(x[i].begin(), x[i].end());
            for (auto p = x[i].begin(); p != t; ++p) {
                int id = -1;
                if (code.count(*p) == 0) {
                    id = code.size();
                    code[*p] = id;
                }
                else
                    id = code[*p];
                if (widx.size() < id + 1) {
                    //cout << *p << " " << id << " " << widx.size() << endl;
                    widx.emplace_back(1, i);
                    word.push_back(*p);
                }
                else
                    widx[id].push_back(i);
            }
            dl[i] = x[i].size();
        }
        int m = widx.size();
        IntegerVector ptr(m + 1); ptr[0] = 0;
        for (int i = 0; i < m; ++i)
            ptr[i + 1] = ptr[i] + widx[i].size();
        int nnz = ptr[m];
        IntegerVector idx(nnz);
        DoubleVector val(nnz, 1.0);
        IntegerVector df(nnz);
        for (int i = 0; i < m; ++i) {
            memcpy(idx.begin() + ptr[i], widx[i].data(), sizeof(int) * widx[i].size());
            df[i] = widx[i].size();
        }
        S4 ans("dtm");
        ans.slot("p") = ptr;
        ans.slot("i") = idx;
        ans.slot("x") = val;
        ans.slot("Dim") = IntegerVector::create(n, m);
        ans.slot("Dimnames") = List::create(R_NilValue, word);
        ans.slot("dl") = dl;
        ans.slot("df") = df;
        return ans;
    } else {
        for (int i = 0; i < n; ++i) {
            for (auto p = x[i].begin(); p != x[i].end(); ++p) {
                int id = -1;
                if (code.count(*p) == 0) {
                    id = code.size();
                    code[*p] = id;
                }
                else
                    id = code[*p];
                if (widx.size() < id + 1) {
                    widx.emplace_back(1, i);
                    word.push_back(*p);
                }
                else
                    widx[id].push_back(i);
            }
            dl[i] = x[i].size();
        }
        int m = widx.size();
        IntegerVector ptr(m + 1); ptr[0] = 0;
        IntegerVector df(m);
        for (int i = 0; i < m; ++i) {
            int uniques = 1;
            for (int j = 1; j < widx[i].size(); ++j)
                if (widx[i][j] != widx[i][j - 1])
                    ++uniques;
            ptr[i + 1] = ptr[i] + uniques;
            df[i] = uniques;
        }
        int nnz = ptr[m];
        IntegerVector idx(nnz);
        DoubleVector val(nnz);
        for (int i = 0; i < m; ++i) {
            int p = ptr[i];
            idx[p] = widx[i][0];
            val[p] = 1;
            for (int j = 1; j < widx[i].size(); ++j) {
                if (widx[i][j] != widx[i][j - 1]) {
                    ++p;
                    idx[p] = widx[i][j];
                    val[p] = 1;
                    if (idf)
                        val[p - 1] *= log(1 + (double)n / df[i]);
                }
                else
                    ++val[p];
            }
            if (idf)
                val[p] *= log(1 + (double)n / df[i]);
        }
        S4 ans("dtm");
        ans.slot("p") = ptr;
        ans.slot("i") = idx;
        ans.slot("x") = val;
        ans.slot("Dim") = IntegerVector::create(n, m);
        ans.slot("Dimnames") = List::create(R_NilValue, word);
        ans.slot("dl") = dl;
        ans.slot("df") = df;
        return ans;
    }
}

RcppExport SEXP C_dtm(SEXP xxSEXP, SEXP ngramSEXP, SEXP tfSEXP, SEXP idfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::vector<std::vector<std::string> >& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< std::vector<int>& >::type ngram(ngramSEXP);
    Rcpp::traits::input_parameter< bool >::type tf(tfSEXP);
    Rcpp::traits::input_parameter< bool >::type idf(idfSEXP);
    rcpp_result_gen = Rcpp::wrap(_C_dtm(xx, ngram, tf, idf));
    return rcpp_result_gen;
END_RCPP
}

void _C_update_dtm_subset(S4 x) {
    IntegerVector ptr = x.slot("p");
    IntegerVector idx = x.slot("i");
    IntegerVector dim = x.slot("Dim");
    int nrow = dim[0], ncol = dim[1];
    IntegerVector df(ncol);
    IntegerVector dl(nrow, 0);
    for (int i = 0; i < ncol; ++i) {
        df[i] = ptr[i + 1] - ptr[i];
        for (int j = ptr[i]; j < ptr[i + 1]; ++j)
            ++dl[idx[j]];
    }
    x.slot("dl") = dl;
    x.slot("df") = df;
}

RcppExport SEXP C_update_dtm_subset(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< S4 >::type x(xSEXP);
    _C_update_dtm_subset(x);
    return R_NilValue;
END_RCPP
}
