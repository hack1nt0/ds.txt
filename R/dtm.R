dtm <- function(x, ngram = 1L:1L, tf = T, idf = T, test) {
    stopifnot(is.list(x) && all(sapply(x, is.character)))
    ngram = as.integer(ngram)
    stopifnot(length(ngram) > 0 && identical(is.unsorted(ngram), FALSE)
              && ngram[[1L]] <= ngram[[length(ngram)]])
    gfrom = ngram[[1L]]
    gto = ifelse(length(ngram) == 1L, ngram[[1L]] + 1, ngram[[length(ngram)]] + 1)
    tf = as.logical(tf)
    idf = as.logical(idf)
    if (missing(test))
        test = NULL
    else
        stopifnot(is.list(test) && all(sapply(test, is.character)))
    ans = .Call('C_dtm', PACKAGE = 'ds.txt', x, gfrom, gto, tf, idf, test)
}

cv.dtm <- function(x, nfolds = 10L, foldid, parallel = FALSE, ...) {
    stopifnot(is.list(x) && all(sapply(x, is.character)))
    nfolds = as.integer(nfolds)
    n = length(x)
    if (missing(foldid)) {
        if (nfolds < 3L || n <= nfolds) {
            warning(paste0("nfolds should be of bound [3, ", n, ")."))
            nfolds = min(n, max(3L, nfolds))
        }
        foldid = sample(rep(seq(nfolds), length = n))
    }
    else {
        foldid = as.integer(foldid)
        stopifnot(length(foldid) == n && max(diff(sort(foldid))) == 1L)
        nfolds = max(foldid)
    }
    ans = as.list(seq(nfolds))
    if (parallel) {
        ans = foreach(i == seq(nfolds), .packages = c("ds.txt")) %dopar% {
            train = foldid == i
            test = !train
            dtm(x[train], ..., test = x[test])
        }
    }
    else {
        for(i in seq(nfolds)) {
            train = foldid == i
            test = !train
            ans[[i]] = dtm(x[train], ..., test = x[test])
        }
    }
    class(ans) = "cv.dtm"
    ans
}

is.dtm <- function(x) class(x) == "dtm"

apply.dtm <- function(x, MARGIN, FUN, ..., simplify = TRUE) {
    stopifnot(inherits(x, "dgCMatrix"))
    FUN = match.fun(FUN)
    MARGIN = as.integer(MARGIN)
    stopifnot((MARGIN == 1L || MARGIN == 2L))
    if (MARGIN == 1L) {
        n = nrow(x)
        ans = vector("list", n)
        for (i in 1:n) {
            ans[[i]] = do.call(FUN, list(x = x[i, , drop = F], ...))
            # if (i %% 10 == 0) print(paste(i, ans[[i]]))
        }
    }
    else {
        n = ncol(x)
        ans = vector("list", n)
        for (i in 1:n)
            ans[[i]] = do.call(FUN, list(x = x[ , i, drop = F], ...))
    }
    if (!identical(simplify, FALSE) && length(ans))
        simplify2array(ans, higher = (simplify == "array"))
    else
        ans
}
