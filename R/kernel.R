crosskern <- function(X, Y = NULL, method = "rbf") {
    methods = c("rbf", "linear", "guassian")
    method = match.arg(method, methods)
    id.method = match(method, methods)
    if (is.null(Y)) Y = X
    stopifnot((is.matrix(X) || inherits(X, "Matrix"))
                  && (is.matrix(Y) || inherits(Y, "Matrix")))
    # .Call("C_crosskern", PACKAGE = "ds.txt", X, Y, id.method)
}
