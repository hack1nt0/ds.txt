fxy <- function(X, Y, potentials = list (
    list(word = dtm(ngram = 1:1, weight = 0L), y = NA, value = 1)
), lower.f = 1L, lower.d = 1L) {
    Y = as.factor(Y)
    levels.Y = levels(Y)
    Y = as.integer(Y)
    np = length(potentials)
    xs = list()
    ws = list()
    for (p in potentials) {
        predicate = p[[1L]]
        y = p[[2L]]
        if (is.na(y))
            y = 1:length(levels.Y)
        else {
            y = as.integer(factor(p[[2L]], levels = levels.Y))
            if (is.na(y))
                stop(paste("Value of feature Function outcome [", P[[2L]], "] not in given levels [",
                           paste(levels.Y, sep = ", "), "]."))
        }
        value = ifelse(is.null(p$value), 1, p$value)
        SM = as(predicate(X), "dgCMatrix") #todo
        prefix = names(p)[[1L]]
        colnames(SM) = str_c(prefix, colnames(SM), sep = "_")
        if (!is.null(xs$prefix))
            stop(paste("Duplicated Feature Functions(name): ", prefix))
        xs$prefix = SM
        W = matrix(0, nrow = ncol(SM), ncol = length(levels.Y))
        W[, y] = 1
        ws$prefix = W
    }
    ans = list()
    FX = do.call(cbind, xs)
    FW = do.call(cbind, ws)
    validWi = (lower.f <= rowSums(FW))
    ans$features = c(sum(validWi), nrow(FW), lower.f)
    FW = FW[validWi,, drop = F]
    FX = FX[, validWi, drop = F]
    validXi = (lower.d <= Matrix::rowSums(FX))
    FX = FX[validXi,, drop = F]
    ans$Y = Y[validXi]
    ans$inputs = c(sum(validXi), length(X), lower.d)
    ans$X = FX
    ans$W = FW
    ans$levels.Y = levels.Y
    ans$potentials = potentials
    ans$call = match.call()
    class(ans) = "fxy"
    ans
}

print.fxy <- function(obj) {
    cat("Call: ")
    print(obj$call)
    stat = matrix(c(obj$features, obj$inputs), nrow = 2L, byrow = T)
    rownames(stat) = c("Feature", "Input")
    colnames(stat) = c("Valid", "Total", "Cutoff(lower bound)")
    print(stat)
}
