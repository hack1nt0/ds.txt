fxy <- function(X, Y, potentials = list (
    list(word = dtm(ngram = 1:1, weight = 0L), y = NA, score = 1)
), lower.f = 1L, lower.d = 1L) {

    stopifnot(is.vector(X) && sum(sapply(X, is.character)) == length(X))
    stopifnot(is.atomic(Y))
    stopifnot(length(X) == length(Y))
    # todo check NA of X and Y
    lower.f = as.integer(lower.f)
    lower.d = as.integer(lower.d)
    class.Y = class(Y)
    Y = as.factor(Y)
    levels.Y = levels(Y)
    Y = as.integer(Y)
    np = length(potentials)
    xs = list()
    ss = list()
    np = length(potentials)
    predicates = vector("list", np)
    ys = vector("numeric", np)
    scores = vector("numeric", np)
    ips = 1:np
    for (ip in ips) {
        p = potentials[[ip]]
        predicate = predicates[[ip]] = p[[1L]]
        y = ys[ip] = p[[2L]]
        if (is.na(y))
            y = 1:length(levels.Y)
        else {
            stopifnot(class(y) == class.Y)
            y = as.integer(factor(p[[2L]], levels = levels.Y))
            if (is.na(y))
                stop(paste0("y of potential functions contains class-tag not belongs to Y, i.e, ", P[[2L]], " not in [", paste(levels.Y, sep = ", "), "]."))
        }
        SM = as(predicate(X), "dgCMatrix") #todo
        prefix = names(p)[[1L]]
        colnames(SM) = str_c(prefix, colnames(SM), sep = "_")
        if (!is.null(xs$prefix))
            stop(paste("Duplicated Feature Functions(name): ", prefix))
        xs$prefix = SM
        S = matrix(0, nrow = ncol(SM), ncol = length(levels.Y))
        score = scores[ip] = ifelse(is.null(p$score), 1, p$score) #todo how to use
        S[, y] = score
        ss$prefix = S
    }
    ans = list()
    FX = do.call(cbind, xs)
    FS = do.call(cbind, ss)
    validWi = (lower.f <= Matrix::colSums(FX))
    ans$features = c(sum(validWi), nrow(FS), lower.f)
    FS = FS[validWi,, drop = F]
    FX = FX[, validWi, drop = F]
    validXi = (lower.d <= Matrix::rowSums(FX))
    FX = FX[validXi,, drop = F]
    ans$Y = Y[validXi]
    ans$inputs = c(sum(validXi), length(X), lower.d)
    ans$FX = FX
    ans$FS = FS
    ans$levels.Y = levels.Y
    ans$potentials = data.table(predicate = predicates, y = ys, score = scores)
    nx = nrow(FX)
    np = ncol(FX)
    stopifnot(ncol(FX) == nrow(FS))
    nc = ncol(FS)
    nf = np * nc
    ans$nx = nx; ans$np = np; ans$nc = nc; ans$nf = nf;
    iNozeros = data.table(Matrix::which(FX != 0, arr.ind = T))
    adjDfm = iNozeros[, list(adj = list(col)), keyby = row]$adj
    adjFdm = iNozeros[, list(adj = list(row)), keyby = col]$adj
    actualMean = matrix(0, nrow = np, ncol = nc)
    for (ip in 1:np) {
        ys = Y[adjFdm[[ip]]]
        for (y in ys)
            actualMean[ip, y] = actualMean[ip, y] + 1
    }
    ans$adjDfm = adjDfm
    ans$adjFdm = adjFdm
    ans$actualMean = actualMean
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
