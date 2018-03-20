#include <Rcpp.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
using namespace Rcpp;
using namespace std;

typedef std::vector<std::vector<std::string> > Text;

RObject _C_dtm(Text& train, int gfrom, int gto,
          bool tf, bool idf, Text& test) {

    List ans(2);
    Text* xp[2] = {&train, &test};
    int n[2] = {(int)train.size(), (int)test.size()};
    int m = 0;
    map<string, int> code;
    vector<string> word;
    IntegerVector* df = NULL;
    for (int t = 0; t < 2; ++t) {
        if (n[t] == 0) continue;
        if (!(gfrom == 1 && gto == 2)) {
            xp[t] = new Text(n[t], vector<string>());
            Text& newx = *xp[t];
            Text& oldx = t == 0 ? train : test;
            for (int i = 0; i < n[t]; ++i) {
                for (int j = 0; j + gfrom <= oldx[i].size(); ++j) {
                    stringbuf sb;
                    for (int k = j; k < j + gfrom; ++k)
                        sb.sputn(oldx[i][k].c_str(), oldx[i][k].size());
                    newx[i].push_back(sb.str()); //todo
                    for (int len = gfrom + 1; len < gto; ++len) {
                        if (j + len > oldx[i].size()) break;
                        sb.sputn(oldx[i][j + len - 1].c_str(), oldx[i][j + len - 1].size());
                        newx[i].push_back(sb.str()); //todo
                    }
                }
            }
        }
        Text& x = *xp[t];
        vector<vector<int> > widx(t == 0 ? 0 : m, vector<int>());
        IntegerVector dl(n[t]);
        if (!tf) {
            for (int i = 0; i < n[t]; ++i) {
                sort(x[i].begin(), x[i].end());
                auto bow = std::unique(x[i].begin(), x[i].end());
                for (auto p = x[i].begin(); p != bow; ++p) {
                    if (t == 0) {
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
                        ++dl[i];
                    }
                    else {
                        if (code.count(*p) > 0) {
                            widx[code[*p]].push_back(i);
                            ++dl[i];
                        }
                    }
                }
            }
            if (t == 0) {
                m = widx.size();
                df = new IntegerVector(m);
            }
            IntegerVector ptr(m + 1); ptr[0] = 0;
            for (int i = 0; i < m; ++i)
                ptr[i + 1] = ptr[i] + widx[i].size();
            int nnz = ptr[m];
            IntegerVector idx(nnz);
            DoubleVector val(nnz, 1.0);
            for (int i = 0; i < m; ++i) {
                memcpy(idx.begin() + ptr[i], widx[i].data(), sizeof(int) * widx[i].size());
                df->operator[](i) = widx[i].size();
            }
            S4 res("dtm");
            res.slot("p") = ptr;
            res.slot("i") = idx;
            res.slot("x") = val;
            res.slot("Dim") = IntegerVector::create(n[t], m);
            res.slot("Dimnames") = List::create(R_NilValue, word);
            res.slot("dl") = dl;
            res.slot("df") = *df;
            ans[t] = res;
        } else {
            for (int i = 0; i < n[t]; ++i) {
                for (auto p = x[i].begin(); p != x[i].end(); ++p) {
                    if (t == 0) {
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
                    else {
                        if (code.count(*p) > 0) {
                            widx[code[*p]].push_back(i);
                            ++dl[i];
                        }
                    }
                }
                dl[i] = x[i].size();
            }
            if (t == 0) {
                m = widx.size();
                df = new IntegerVector(m);
            }
            IntegerVector ptr(m + 1); ptr[0] = 0;
            for (int i = 0; i < m; ++i) {
                int uniques = 1;
                for (int j = 1; j < widx[i].size(); ++j)
                    if (widx[i][j] != widx[i][j - 1])
                        ++uniques;
                    ptr[i + 1] = ptr[i] + uniques;
                    df->operator[](i) = uniques;
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
                            val[p - 1] *= log(1 + (double)n[t] / df->operator[](i));
                    }
                    else
                        ++val[p];
                }
                if (idf)
                    val[p] *= log(1 + (double)n[t] / df->operator[](i));
            }
            S4 res("dtm");
            res.slot("p") = ptr;
            res.slot("i") = idx;
            res.slot("x") = val;
            res.slot("Dim") = IntegerVector::create(n[t], m);
            res.slot("Dimnames") = List::create(R_NilValue, word);
            res.slot("dl") = dl;
            res.slot("df") = *df;
            ans[t] = res;
        }
    }
    return test.size() == 0 ? ans[0] : ans;
}

RcppExport SEXP C_dtm(SEXP xxSEXP, SEXP ngramSEXP, SEXP gfromSEXP, SEXP gtoSEXP, SEXP tfSEXP, SEXP idfSEXP
                          , SEXP testSEXP
) {
    Rcpp::RNGScope rcpp_rngScope_gen;
    BEGIN_RCPP
        Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< std::vector<std::vector<std::string> >& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< int >::type gfrom(gfromSEXP);
    Rcpp::traits::input_parameter< int >::type gto(gtoSEXP);
    Rcpp::traits::input_parameter< bool >::type tf(tfSEXP);
    Rcpp::traits::input_parameter< bool >::type idf(idfSEXP);
    Text test = Rf_isNull(testSEXP) ? Text() : as<Text>(testSEXP);
    rcpp_result_gen = Rcpp::wrap(_C_dtm(xx, gfrom, gto, tf, idf, test));
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
