
fxy <- function(X, Y = NULL, ngram = 1L, tf = F, idf = F, ff = list(
  x.equal.y = function(i) identical(X[[i]], Y[[i]])
  , row.idx = function(i) i
  ))
{
  stopifnot(is.list(X) && all(sapply(X, is.character)))
  require(Matrix)
  if (ngram <= 0L)
    A <- NULL
  else {
    A <- dtm(X, tf = as.logical(tf), idf = as.logical(idf))
    warning("X is not list of character', entering test mode...")
    A <- rsparsematrix(length(X), 2, 0.5)
    colnames(A) = c("T1", "T2")
  }
  B0 = vector("list", length(ff))
  names.B0 = vector("character", length(ff))
  j <- 1L
  names.ff = names(ff)
  for (fi in 1:length(ff)) {
    f = ff[[fi]]
    col = vector("character", length(X))
    for (i in 1:length(X)) col[[i]] = f(i)
    col = as.factor(col)
    if (nlevels(col) >= 2L) {
      B0[[j]] = col
      names.B0[[j]] = ifelse(is.null(names.ff[[fi]]),
                             as.character(enquote(body(f)))[[2L]],
                             names.ff[[fi]])
      j = j + 1L
    }
    else
      warning(paste0("Feature function(", names.ff[[fi]], " = ", as.character(enquote(body(f)))[[2L]],
                     ")must be meaningful(e.g. nlevels(extracted.features) >= 2L).\n"))
  }
  if (j > 1L) {
    names(B0) = names.B0
    B <- sparse.model.matrix(~., as.data.frame(B0[1:(j - 1L)]))
    cbind(B, A)
  }
  else {
    C = cbind(1, A)
    colnames(C)[[1L]] = "(Intercept)"
    C
  }
}

#fxy(list("A", "B"))

dtm <- function(x, ngram = 1L:1L, tf = T, idf = T) {
    stopifnot(length(ngram) <= 2 && is.list(x) && all(sapply(x, is.character)))
    #C_dtm(x, sort(as.integer(ngram)), tf, idf)
    .Call('ds_txt_C_dtm', PACKAGE = 'ds.txt', x, as.integer(ngram), tf, idf)
}
