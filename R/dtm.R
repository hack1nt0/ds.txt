dtm <- function(x, ngram = 1L:1L, tf = T, idf = T) {
    stopifnot(is.list(x) && all(sapply(x, is.character)))
    ngram = as.integer(ngram)
    stopifnot(length(ngram) > 0 && identical(is.unsorted(ngram), FALSE))
    ngram = ifelse(length(ngram) == 1L, c(ngram[[1L]], ngram[[1L]] + 1L),
                   c(ngram[[1L]], ngram[[length(ngram)]] + 1L))
    ans = .Call('C_dtm', PACKAGE = 'ds.txt', x, ngram, tf, idf)
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
