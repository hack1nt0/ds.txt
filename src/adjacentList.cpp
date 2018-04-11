#include <Rcpp.h>
#include <vector>
#include <map>
using namespace Rcpp;
using namespace std;

List _C_new_adjacent_list_dgCMatrix(int nrow, IntegerVector& ptr, IntegerVector& idx, CharacterVector& colnames, CharacterVector& dict) {
    int ncol = ptr.size() - 1;
    map<String, int> word2idx;
    int id = 1;
    for (String word : dict) {
        if (word2idx.count(word) > 0)
            Rcpp::stop("Duplicated word in dict.");
        word2idx[word] = id++;
    }
    vector<vector<int> > adj(nrow, vector<int>());
    for (int c = 0; c < ncol; ++c) {
        if (word2idx.count(colnames[c]) == 0)
            continue;
        int cId = word2idx[colnames[c]];
        for (int p = ptr[c]; p < ptr[c + 1]; ++p) {
            int r = idx[p];
            adj[r].push_back(cId);
        }
    }
    return wrap(adj);
}

RcppExport SEXP C_new_adjacent_list_dgCMatrix(SEXP XSEXP, SEXP dictSEXP) {
    Rcpp::RNGScope rcpp_rngScope_gen;
    BEGIN_RCPP
        Rcpp::RObject rcpp_result_gen;
    S4 X(XSEXP);
    int nrow = as<IntegerVector>(X.slot("Dim"))[0];
    IntegerVector ptr = X.slot("p");
    IntegerVector idx = X.slot("i");
    CharacterVector colnames = as<List>(X.slot("Dimnames"))[1];
    CharacterVector dict = as<CharacterVector>(dictSEXP);
    rcpp_result_gen = Rcpp::wrap(_C_new_adjacent_list_dgCMatrix(nrow, ptr, idx, colnames, dict));
    return rcpp_result_gen;
    END_RCPP
}

List _C_new_adjacent_list_integer(IntegerVector& X, int ngroup) {
    vector<vector<int> > adj(ngroup, vector<int>());
    int nx = X.size();
    for (int i = 0; i < nx; ++i)
        adj[X[i] - 1].push_back(i + 1);
    return wrap(adj);
}

RcppExport SEXP C_new_adjacent_list_integer(SEXP XSEXP, SEXP ngroupSEXP) {
    Rcpp::RNGScope rcpp_rngScope_gen;
    BEGIN_RCPP
        Rcpp::RObject rcpp_result_gen;
    IntegerVector X(XSEXP);
    int ngroup = as<int>(ngroupSEXP);
    rcpp_result_gen = Rcpp::wrap(_C_new_adjacent_list_integer(X, ngroup));
    return rcpp_result_gen;
    END_RCPP
}
