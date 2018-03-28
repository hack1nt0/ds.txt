me <- function(X, Y, niters = 100L, verbose = TRUE,
               predicates = list(
                   function(X, Y, i) str_detect(X[[i]], "澳门"),
                   function(X, Y, i) str_detect(X[[i]], "余额"),
                   function(X, Y, i) str_detect(X[[i]], "OX")
               ),
               y = c(TRUE, FALSE, FALSE),
               lower = 1L, upper = .Machine$integer.max
) {
    stopifnot(is.vector(X) && sum(sapply(X, is.character)) == length(X))
    stopifnot(class(Y) == class(y) && is.atomic(Y))
    stopifnot(length(X) == length(Y))
    stopifnot(length(predicates) == length(y))
    if (sum(sapply(X, is.null)) > 0L || sum(sapply(Y, is.null)) > 0L
        || sum(sapply(X, is.na)) > 0L || sum(sapply(Y, is.na)) > 0L) {
        stop("Input data contains either NA or NULL.")
    }
    if (sum(sapply(y, is.null)) > 0L || sum(sapply(y, is.null)) > 0L)
        stop("y of feature functions contain either NA or NULL. You can use file input to illustrate\n
             class-unrelated predicates.")
    Y = as.factor(Y)
    levels.Y = levels(Y)
    Y = as.integer(Y)
    y = as.integer(factor(y, levels = levels.Y))
    if (sum(sapply(y, is.na)) > 0L)
        stop("y of potential functions contains class-tag not belongs to Y.")
    if (0 < lower && lower < 1) lower = as.integer(n * lower)
    if (0 < upper && upper < 1) upper = as.integer(n * upper)

    n = length(X); m = length(levels.Y)
    np = length(predicates)
    # weight = vector("numeric", np)
    weight = rnorm(np) # todo
    P = matrix(0, nrow = n, ncol = m) # probability
    M = vector("list", n)
    actualMean = rep(0, np)
    maxentMean = rep(0, np)
    whichP = rep(TRUE, np)
    ips = 1:np
    ixs = 1:n
    for (ix in ixs) {
        for (ip in ips) {
            whichP[ip] = predicates[[ip]](X, Y, ix)
            if (whichP[ip] && y[ip] == Y[ix])
                actualMean[ip] = actualMean[ip] + 1
        }
        M[[ix]] = ips[whichP]
    }
    valid <- function(ip) {
        lower <= actualMean[ip] & actualMean[ip] < upper
    }
    for (ix in ixs) {
        M[[ix]] = M[[ix]][valid(M[[ix]])]
    }
    ips = ips[valid(ips)]
    losses = rep(NA, niters)
    it = 1L
    loss <- function(weight) {
        lhs = 0; rhs = 0
        for (ix in ixs) {
            for (ip in M[[ix]]) {
                yip = y[ip]
                P[ix, yip] = P[ix, yip] + weight[yip]
                if (yip == Y[ix])
                    rhs = rhs + weight[yip]
            }
        }
        lhs = log(rowSums(exp(P))) #todo
        ans = sum(lhs) - rhs
        losses[it] = ans
        it = it + 1L
        # print(paste("loss: ", ans))
        ans
    }
    grad <- function(weight) {
        for (ix in ixs) {
            for (ip in M[[ix]]) {
                yip = y[ip]
                P[ix, yip] = P[ix, yip] + weight[yip]
            }
        }
        P[,] = exp(P)
        P[,] = P / rowSums(P) #todo
        maxentMean[] = 0
        for (ix in ixs) {
            for (ip in ips) {
                maxentMean[ip] = maxentMean[ip] + P[ix, y[ip]]
            }
        }
        ans = maxentMean - actualMean
        # print(paste("grad: ", ans))
        ans
    }
    opt.obj = optim(weight, fn = loss, gr = grad, method = "L-BFGS", control = list(
        trace = ifelse(verbose, 10L, 0L), maxit = niters))
    ans = list(opt.obj = opt.obj, predicates = predicates, y = y,
               ips = ips, M = M, losses = losses,
               levels.Y = levels.Y, call = match.call())
    class(ans) = "me"
    ans
}

predict.me <- function(obj, X, ...) {
    nx = length(X)
    nclass = length(obj$levels.Y)
    P = matrix(0, nx, nclass)
    maxentY = rep(0, nx) # class-id start from 1L
    predicates = obj$predicates
    y = obj$y
    weight = obj$opt.obj$par
    ips = 1:length(predicates)
    for (ix in 1:nx) {
        for (ip in ips)
            if (predicates[[ip]](X, maxentY, ix))
                P[ix, y[ip]] = P[ix, y[ip]] + weight[ip]
    }
    P[,] = exp(P)
    P[,] = P / rowSums(P)
    ans = list(P = P, Y = obj$levels.Y[apply(P, 1L, which.max)])
    class(ans) = "predict.me"
    ans
}

plot.me <- function(obj, ...) {
    plot(obj$losses, xlab = "Iteration", ylab = "Negative Log Likelihood")
}
