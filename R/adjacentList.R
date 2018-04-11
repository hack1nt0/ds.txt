adjacent.list <- function(X, ...) UseMethod("adjacent.list")

adjacent.list.dgCMatrix <- function(X, byrow = T, dict = NULL) {
    if (is.null(dict)) {
        iNozeros = data.table(Matrix::which(X != 0, arr.ind = T))
        ans =
            if (byrow)
                iNozeros[, list(adj = list(col)), keyby = row]$adj
            else
                iNozeros[, list(adj = list(row)), keyby = col]$adj
    } else {
        stopifnot(is.character(dict))
        if (!byrow)
            X = t(X)
        ans = .Call("C_new_adjacent_list_dgCMatrix", PACKAGE = "ds.txt", X, dict)
    }
    class(ans) = "adjacent.list"
    ans
}

adjacent.list.integer <- function(X, byrow = T, ngroup = max(X)) {
    ans = NULL
    if (byrow)
        ans = as.list(X)
    else {
        stopifnot(max(X) < ngroup)
        ans = .Call("C_new_adjacent_list_integer", PACKAGE = "ds.txt", X, ngroup)
    }
    class(ans) = "adjacent.list"
    ans
}
