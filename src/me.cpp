#include <Rcpp.h>
#include <vector>
#include <string>
#include <map>
#include <algorithm>
#include <cstring>
#include <unordered_map>
#include <iostream>
using namespace Rcpp;
using namespace std;

template <class V>
double logsumexp(const V & v) {
    double lse = 0;
    double maxV = v[0];
    for (auto p = v.begin(); p != v.end(); ++p)
        maxV = max(maxV, *p);
    for (auto p = v.begin(); p != v.end(); ++p)
        lse += exp(*p - maxV);
    lse = maxV + log(lse);
    return lse;
}

NumericMatrix makeFC(int np, int nc, IntegerVector  ptrP, IntegerVector  indX, IntegerVector Y) {
    NumericMatrix FC(np, nc); FC.fill(.0);
    for (int ip = 0; ip < np; ++ip) {
        for (int pp = ptrP[ip]; pp < ptrP[ip + 1]; ++pp) {
            int y = Y[indX[pp]] - 1;
            FC(ip, y) += 1;
        }
    }
    return FC;
}

SEXP _C_me(int nx, int np, int nc,
        IntegerVector  ptrX, IntegerVector  indP,
        IntegerVector  ptrP, IntegerVector  indX,
        IntegerVector  Y, double lambda, double alpha,
        NumericMatrix  FW, NumericMatrix  FS,
        int niter, double tolerance, double eta, double beta1, double beta2, double epsilon, bool AMSGrad,
        bool verbose,
        NumericMatrix  FC, NumericMatrix  P,
        NumericMatrix  meanGrad,
        NumericMatrix  varGrad,
        NumericMatrix  maxVarGrad,
        bool keep
        ) {
    cout << "eta: " << eta << endl;
    cout << beta1 << endl;
    cout << beta2 << endl;
    meanGrad.fill(0.);
    varGrad.fill(0.);
    maxVarGrad.fill(0.);
    vector<double> losses;
    double powBeta1 = beta1;
    double powBeta2 = beta2;
    LogicalMatrix converaged(np, nc); //todo
    converaged.fill(false);
    int nconveraged = 0, nf = np * nc;
    for (int it = 0, ib = 0; it < niter; ++it) {
        if (nconveraged == nf)
            break;
        // loss
        P.fill(0.);
        double lhs, rhs;
        lhs = rhs = 0.;
        for (int ix = 0; ix < nx; ++ix) {
            for (int px = ptrX[ix]; px < ptrX[ix + 1]; ++px) {
                int ip = indP[px];
                for (int y = 0; y < nc; ++y)
                    P(ix, y) += FW(ip, y) * FS(ip, y);
            }
            rhs += P(ix, Y[ix] - 1);
            double lse = logsumexp(P.row(ix));
            lhs += lse;
            for (int y = 0; y < nc; ++y)
                P(ix, y) = exp(P(ix, y) - lse);
        }
        double loss = (lhs - rhs) / nx;
        losses.push_back(loss);
        if (verbose)
            Rcout << "Epoch " << it << "\t, loss " << loss << endl;
        if (losses.size() >= 2 && abs(losses[it] - losses[it - 1]) <= tolerance)
            break;

        //grad
        for (int ip = 0; ip < np; ++ip) {
            for (int y = 0; y < nc; ++y) {
                if (converaged(ip, y))
                    continue;
                double score = FS(ip, y);
                if (score <= 0.)
                    continue;
                double fc = 0;
                for (int pp = ptrP[ip]; pp < ptrP[ip + 1]; ++pp) {
                    int ix = indX[pp];
                    fc += P(ix, y);
                }
                double oldw = FW(ip, y);
                double signOldW = oldw < -tolerance ? -1 : (oldw > +tolerance ? +1 : 0);
                double g = (fc - FC(ip, y)) * score / nx + lambda * ((1 - alpha) * 2. * oldw + alpha * signOldW);
                // double g = (fc - FC(ip, y)) * score / nx;
                double mean = beta1 * meanGrad(ip, y) + (1. - beta1) * g;
                double var  = beta2 * varGrad(ip, y) + (1. - beta2) * g * g;
                double maxVar = max(maxVarGrad(ip, y), var);
                double neww = 0;
                if (AMSGrad)
                    neww = oldw - eta * (mean) / (sqrt(maxVar) + epsilon);
                else
                    neww = oldw - eta * (mean / (1 - powBeta1)) / (sqrt(var / (1 - powBeta2)) + epsilon);
                if (abs(neww - oldw) < tolerance) {
                    converaged(ip, y) = true;
                    ++nconveraged;
                    continue;
                }
                FW(ip, y) = neww;
                meanGrad(ip, y) = mean;
                varGrad(ip, y) = var;
                maxVarGrad(ip, y) = maxVar;
            }
        }
        if (it > 0 && ib == 0) {
            powBeta1 *= beta1;
            powBeta2 *= beta2;
        }
    }
    if (keep) {
        return List::create(Named("losses") = losses,
                            Named("mean.gr") = meanGrad,
                            Named("var.gr") = varGrad,
                            Named("max.var.gr") = maxVarGrad,
                            Named("P") = P,
                            Named("converaged") = converaged
                            );
    } else
        return wrap(losses);
}


SEXP _C_cv_me(
        IntegerVector ptrTrainX, IntegerVector indTrainP,
        IntegerVector ptrTrainP, IntegerVector indTrainX,
        IntegerVector trainY,
        IntegerVector ptrTestX, IntegerVector indTestP,
        IntegerVector testY,
        NumericVector lambdas, NumericVector alphas,
        String measureType,
        NumericMatrix initFW, NumericMatrix FS, int niter, double tolerance,
        double eta, double beta1, double beta2, double epsilon, bool AMSGrad,
        bool verbose) {
    int nTrainX = trainY.size();
    int nTestX = ptrTestX.size() - 1;
    int np = initFW.rows();
    int nc = initFW.cols();
    int nlambda = lambdas.size();
    int nalpha = alphas.size();
    vector<double> losses;
    NumericMatrix meanGrad(np, nc);
    NumericMatrix varGrad(np, nc);
    NumericMatrix maxVarGrad(np, nc);
    NumericMatrix trainP(nTrainX, nc);
    NumericMatrix testP(nTestX, nc);
    NumericMatrix FC = makeFC(np, nc, ptrTrainP, indTrainX, trainY);
    NumericMatrix FW; //todo
    for (int ia = 0; ia < nalpha; ++ia) {
        double alpha = alphas[ia];
        FW = initFW;
        for (int il = 0; il < nlambda; ++il) {
            double lambda = lambdas[il];
            _C_me(nTrainX, np, nc, ptrTrainX, indTrainP, ptrTrainP, indTrainX, trainY, lambda, alpha, FW, FS, niter, tolerance,
                 eta, beta1, beta2, epsilon, AMSGrad, verbose,
                 FC, trainP, meanGrad, varGrad, maxVarGrad, false);
            testP.fill(.0);
            for (int ix = 0; ix < nTestX; ++ix) {
                for (int pp = ptrTestX[ix]; pp < ptrTestX[ix + 1]; ++pp) {
                    int ip = indTestP[pp];
                    for (int ic = 0; ic < nc; ++ic)
                        testP(ix, ic) += FW(ip, ic) * FS(ip, ic);
                }
            }
            double loss = 0.;
            if (measureType == "class") {
                for (int ix = 0; ix < nTestX; ++ix) {
                    int guessY = 0;
                    for (int ic = 1; ic < nc; ++ic)
                        if (testP(ix, ic) > testP(ix, guessY))
                            guessY = ic;
                    loss += guessY != testY[ix];
                }
                loss /= nTestX;
            }
            else if (measureType == "log.likelihood") {
                for (int ix = 0; ix < nTestX; ++ix)
                    loss += testP(ix, testY[ix]) - logsumexp(testP.row(ix));
            }
            else if (measureType == "auc") {
                Rcpp::stop("Type auc not implemented yet.");
            }
            losses.push_back(loss);
        }
    }
    return wrap(losses);
}

// C_me
RcppExport SEXP C_me(SEXP nxSEXP, SEXP npSEXP, SEXP ncSEXP,
                     SEXP ptrXSEXP, SEXP indPSEXP, SEXP ptrPSEXP, SEXP indXSEXP, SEXP YSEXP,
                     SEXP lambdaSEXP, SEXP alphaSEXP, SEXP FWSEXP, SEXP FSSEXP,
                     SEXP niterSEXP, SEXP toleranceSEXP, SEXP etaSEXP, SEXP beta1SEXP, SEXP beta2SEXP, SEXP epsilonSEXP, SEXP AMSGradSEXP, SEXP verboseSEXP, SEXP keepSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nx(nxSEXP);
    Rcpp::traits::input_parameter< int >::type np(npSEXP);
    Rcpp::traits::input_parameter< int >::type nc(ncSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ptrX(ptrXSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indP(indPSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ptrP(ptrPSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indX(indXSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type Y(YSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type FW(FWSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type FS(FSSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< bool >::type AMSGrad(AMSGradSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    Rcpp::traits::input_parameter< bool >::type keep(keepSEXP);
    NumericMatrix FC = makeFC(np, nc, ptrP, indX, Y);
    NumericMatrix P(nx, nc);
    NumericMatrix meanGrad(np, nc);
    NumericMatrix varGrad(np, nc);
    NumericMatrix maxVarGrad(np, nc);
    rcpp_result_gen = Rcpp::wrap(_C_me(nx, np, nc, ptrX, indP, ptrP, indX, Y, lambda, alpha, FW, FS, niter, tolerance, eta, beta1, beta2, epsilon, AMSGrad, verbose, FC, P, meanGrad, varGrad, maxVarGrad, keep));
    return rcpp_result_gen;
END_RCPP
}

// C_cv_me
RcppExport SEXP C_cv_me(SEXP ptrTrainXSEXP, SEXP indTrainPSEXP, SEXP ptrTrainPSEXP, SEXP indTrainXSEXP, SEXP trainYSEXP, SEXP ptrTestXSEXP, SEXP indTestPSEXP, SEXP testYSEXP, SEXP lambdasSEXP, SEXP alphasSEXP, SEXP measureTypeSEXP, SEXP initFWSEXP, SEXP FSSEXP, SEXP niterSEXP, SEXP toleranceSEXP, SEXP etaSEXP, SEXP beta1SEXP, SEXP beta2SEXP, SEXP epsilonSEXP, SEXP AMSGradSEXP, SEXP verboseSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< IntegerVector >::type ptrTrainX(ptrTrainXSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indTrainP(indTrainPSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ptrTrainP(ptrTrainPSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indTrainX(indTrainXSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type trainY(trainYSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type ptrTestX(ptrTestXSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type indTestP(indTestPSEXP);
    Rcpp::traits::input_parameter< IntegerVector >::type testY(testYSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type alphas(alphasSEXP);
    Rcpp::traits::input_parameter< String >::type measureType(measureTypeSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type initFW(initFWSEXP);
    Rcpp::traits::input_parameter< NumericMatrix >::type FS(FSSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type eta(etaSEXP);
    Rcpp::traits::input_parameter< double >::type beta1(beta1SEXP);
    Rcpp::traits::input_parameter< double >::type beta2(beta2SEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< double >::type tolerance(toleranceSEXP);
    Rcpp::traits::input_parameter< bool >::type AMSGrad(AMSGradSEXP);
    Rcpp::traits::input_parameter< bool >::type verbose(verboseSEXP);
    rcpp_result_gen = Rcpp::wrap(_C_cv_me(ptrTrainX, indTrainP, ptrTrainP, indTrainX, trainY, ptrTestX, indTestP, testY, lambdas, alphas, measureType, initFW, FS, niter, tolerance, eta, beta1, beta2, epsilon, AMSGrad, verbose));
    return rcpp_result_gen;
END_RCPP
}
