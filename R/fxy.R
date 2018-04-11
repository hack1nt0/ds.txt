fxy <- function(X, Y, potentials = list (
    list(word = ngrams(n = 1:1, tokenize = TRUE), y = NA, score = 1),
    list(nsep = function(X) str_count(X, "[,.，。!！]"), y = NA, score = 1)
),
lower.f = 1L, lower.d = 1L) {
    stopifnot(all(!is.na(X)) & all(!is.null(X)))
    stopifnot(is.vector(X) && sum(sapply(X, is.character)) == length(X))
    ans = list()
    np = length(potentials)
    stopifnot(is.atomic(Y))
    stopifnot(length(X) == length(Y))
    # todo check NA of X and Y
    lower.f = as.integer(lower.f)
    lower.d = as.integer(lower.d)
    class.Y = class(Y)
    Y = as.factor(Y)
    levels.Y = levels(Y)
    nc = length(levels.Y)
    Y = as.integer(Y)
    name.predicate = vector("character", np)
    predicates = vector("list", np)
    ys = vector("list", np)
    scores = vector("numeric", np)
    for (ip in 1:np) {
        p = potentials[[ip]]
        name = names(p)[[1L]]
        if (length(intersect(name.predicate, name)) > 0)
            stop(paste("Duplicated Feature Functions(name): ", prefix))
        name.predicate[ip] = name
        predicate = predicates[[ip]] = p[[1L]]
        y = p[[2L]]
        if (is.na(y))
            y = 1:length(levels.Y)
        else {
            stopifnot(class(y) == class.Y)
            y = as.integer(factor(p[[2L]], levels = levels.Y))
            if (is.na(y))
                stop(paste0("y of potential functions contains class-tag not belongs to Y, i.e, ", P[[2L]], " not in [", paste(levels.Y, sep = ", "), "]."))
        }
        ys[[ip]] = y
        scores[ip] = ifelse(is.null(p$score), 1, p$score) #todo how to use
    }
    fx.obj = fx(X, predicates = predicates, names = name.predicate)
    FX = fx.obj$fx
    FS = .fs(row.extends = fx.obj$ncols, col.idx = ys, scores = scores, nclass = nc)
    stopifnot(ncol(FX) == nrow(FS))
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
    ans$potentials = data.table(name = name.predicate, predicate = predicates, y = ys, score = scores)
    ans$call = match.call()
    class(ans) = "fxy"
    ans
}

fx <- function(X, predicates, names = NULL) {
    np = length(predicates)
    xs = vector("list", np)
    for (ip in 1:np) {
        predicate = predicates[[ip]]
        SM = predicate(X)
        if (!inherits(SM, "dgCMatrix")) {
            if (!is.vector(SM))
                warning(paste0("Predicator ", ip, " may produced ill-formated data."))
            if (length(SM) != length(X))
                stop("You must apply a multiple-input-enabled predicator.")
            SM = ngrams(n = 1, tokenize = FALSE)(lapply(as.vector(SM), as.character))
        }
        SM = as(SM, "dgCMatrix") #todo
        colnames(SM) = str_c(names[ip], colnames(SM), sep = "_")
        xs[[ip]] = SM
    }
    ans = list()
    ans$fx = do.call(cbind, xs)
    ans$ncols = sapply(xs, ncol)
    class(ans) = "fx"
    ans
}

.fs <- function(row.extends, col.idx, scores, nclass) {
    np = length(scores)
    ss = vector("list", np)
    print(row.extends)
    for (ip in 1:np) {
        S = matrix(0, nrow = row.extends[ip], ncol = nclass)
        S[, col.idx[[ip]]] = scores[ip]
        ss[[ip]] = S
    }
    do.call(rbind, ss)
}

ngrams <- function(n = 1:1, tokenize = FALSE, ...) {
    function(X) {
        if (tokenize) X = tk(X, ...)
        dtm(x = X, ngram = n, weight = 0L)
    }
}

is.fxy <- function(x) {
    inherits(x, "fxy")
}

print.fxy <- function(obj) {
    cat("Call: ")
    print(obj$call)
    stat = matrix(c(obj$features, obj$inputs), nrow = 2L, byrow = T)
    rownames(stat) = c("Feature", "Input")
    colnames(stat) = c("Valid", "Total", "Cutoff(lower bound)")
    print(stat)
}
