#include <Rcpp.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
#include <unordered_map>
using namespace Rcpp;
using namespace std;

typedef std::vector<std::vector<std::string> > Text;

List _C_dtm(Text& train, int gfrom, int gto,
          int weight, Text& test) {

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
            int uniques = 0;
            for (int j = 0; j < widx[i].size(); ++j)
                if (j == 0 || widx[i][j] != widx[i][j - 1])
                    ++uniques;
                ptr[i + 1] = ptr[i] + uniques;
                if (t == 0)
                    (*df)[i] = uniques;
        }
        int nnz = ptr[m];
        IntegerVector idx(nnz);
        DoubleVector val(nnz);
        double* dMaxTf = NULL;
        int* dMaxDf = NULL;
        for (int i = 0; i < m; ++i) {
            if (widx[i].size() == 0) continue;
            int p = ptr[i];
            idx[p] = widx[i][0];
            val[p] = 1;
            for (int j = 1; j < widx[i].size(); ++j) {
                if (widx[i][j] != widx[i][j - 1]) {
                    ++p;
                    idx[p] = widx[i][j];
                    val[p] = 1;
                }
                else ++val[p];
            }
            for (int j = ptr[i]; j < ptr[i + 1]; ++j) {
                double tf = 0;
                switch(weight / 10) {
                    case 0:
                        tf = 1; break;
                    case 1:
                        tf = val[j]; break;
                    case 2:
                        tf = val[j] / dl[idx[j]]; break;
                    case 3:
                        tf = 1.0 + log(val[j]); break;
                    case 4:
                        if (dMaxTf == NULL) {
                            dMaxTf = new double[n[t]];
                            memset(dMaxTf, 0, sizeof(double) * n[t]);
                            for (int k = 0; k < m; ++k)
                                for (int p = ptr[k]; p < ptr[k + 1]; ++p)
                                    dMaxTf[idx[p]] = max(dMaxTf[idx[p]], val[p]);
                        }
                        tf = 0.5 + 0.5 * val[j] / dMaxTf[idx[j]]; break;
                    default:
                        throw std::invalid_argument("TF-IDF weighting schemes unvalid.");
                }
                double idf = (*df)[i];
                switch(weight % 10) {
                    case 0:
                        idf = 1; break;
                    case 1:
                        idf = log(n[0] / (1.0 + idf)); break;
                    case 2:
                        idf = log(1.0 + n[0] / (1.0 + idf)); break;
                    case 3:
                        if (dMaxDf == NULL) {
                            dMaxDf = new int[n[t]];
                            memset(dMaxTf, 0, sizeof(int) * n[t]);
                            for (int k = 0; k < m; ++k)
                                for (int p = ptr[k]; p < ptr[k + 1]; ++p)
                                    dMaxDf[idx[p]] = max(dMaxDf[idx[p]], (*df)[k]);
                        }
                        idf = log(dMaxDf[idx[j]] / (1.0 + idf)); break;
                    case 4:
                        idf = log((n[0] - idf + 1.0) / idf); break;
                    default:
                        throw std::invalid_argument("TF-IDF weighting schemes unvalid.");
                }
                val[j] = tf * idf;
            }
        }
        if (dMaxTf != NULL) delete[] dMaxTf;
        if (dMaxDf != NULL) delete[] dMaxDf;
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
    return ans;
}

RcppExport SEXP C_dtm(SEXP xxSEXP, SEXP gfromSEXP, SEXP gtoSEXP, SEXP weightSEXP, SEXP testSEXP) {
    Rcpp::RNGScope rcpp_rngScope_gen;
    BEGIN_RCPP
        Rcpp::RObject rcpp_result_gen;
    Rcpp::traits::input_parameter< Text& >::type xx(xxSEXP);
    Rcpp::traits::input_parameter< int >::type gfrom(gfromSEXP);
    Rcpp::traits::input_parameter< int >::type gto(gtoSEXP);
    Rcpp::traits::input_parameter< int >::type weight(weightSEXP);
    Text test = Rf_isNull(testSEXP) ? Text() : Rcpp::traits::input_parameter< Text& >::type(testSEXP);
    rcpp_result_gen = Rcpp::wrap(_C_dtm(xx, gfrom, gto, weight, test));
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
