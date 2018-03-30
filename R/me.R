me <- function(X, niters = 100L, verbose = TRUE, R.optim = FALSE, nBatches = 2L) {
    library(data.table)
    stopifnot(class(X) == "fxy")
    Y = X$Y; FX = X$FX; FS = X$FS
    nx = X$nx; nc = X$nc; np = X$np; nf = X$nf
    ixs = 1:nx; ips = 1:np
    adjDfm = X$adjDfm; adjFdm = X$adjFdm
    actualMean = X$actualMean
    model <- new.env(parent = emptyenv())
    P = matrix(0, nrow = nx, ncol = nc, dimnames = list(NULL, X$levels.Y)) # probability
    FW = matrix(rnorm(nf), nrow = np, ncol = nc, dimnames = list(colnames(FX),  X$levels.Y))
    FW[FS == 0] = 0
    model$losses = rep(NA, niters)
    model$it = 1L
    loss <- function(weight) {
        FW = matrix(weight, nrow = np, ncol = nc)
        lhs = 0
        rhs = 0
        for (ix in ixs) {
            ifs = adjDfm[[ix]]
            lse = logsumexp(colSums(FS[ifs,, drop = F] * FW[ifs,, drop = F])) #todo
            lhs = lhs + lse
            rhs = rhs + sum(FS[ifs, Y[ix]] * FW[ifs, Y[ix]])
        }
        curLoss = (lhs - rhs) / nx
        model$losses[model$it] = curLoss
        print(paste(model$it, ":", curLoss))
        model$it = model$it + 1L
        curLoss
    }
    grad <- function(weight) {
        FW = matrix(weight, nrow = np, ncol = nc)
        for (ix in ixs) {
            ifs = adjDfm[[ix]]
            P[ix,] = colSums(FS[ifs,, drop = F] * FW[ifs,, drop = F]) #todo
            P[ix,] = exp(P[ix,] - logsumexp(P[ix,]))
        }
        # P[,] = exp(P)
        # P[,] = P / Matrix::rowSums(P)
        # P[,] = P / rowSums(P) #todo necessary?
        for (ip in ips) {
            ids = adjFdm[[ip]]
            FW[ip,] = colSums(P[ids,, drop = F])
        }
        FW[,] = (FW * FS - actualMean) / nx
        as.double(FW) # todo
    }
    # check.grad(loss, grad, FW, epsilon = 1e-9)
    if (R.optim) {
        opt.obj = optim(as.double(FW), fn = loss, gr = grad, method = "L-BFGS", control = list(
            trace = ifelse(verbose, 10L, 0L), maxit = niters))
        FW[,] = opt.obj$par
    } else {
        # Adam
        lambda = 0.1
        beta.mean = 0.9
        beta.var = 0.999
        epsilon = 1e-8
        mean.grad = rep(0, nf)
        var.grad = rep(0, nf)
        make.batches <- function(x, nBatches = 10L) {
            stopifnot(is.integer(x))
            ans = vector("list", nBatches)
            if (length(x) == 1L)  {

            }
        }
        batches = make.batches(nx, nBatches)
        for (it in 1:niters) {
            if (it == 1 || (it - 1) %% nBatches == 0) {
                loss(as.double(FW))
            }
            for (ixs in batches) {
                for (ix in ixs) {
                    ifs = adjDfm[[ix]]
                    P[ix,] = colSums(FS[ifs,, drop = F] * FW[ifs,, drop = F]) #todo
                    P[ix,] = exp(P[ix,] - logsumexp(P[ix,]))
                }
                cur.grad = grad(as.double(FW))
                mean.grad[] = beta.mean * mean.grad + (1 - beta.mean) * cur.grad
                mean.var[]  = beta.var * mean.var + (1 - beta.var) * cur.grad * cur.grad
                FW[,] = FW - lambda * mean.grad / (sqrt(mean.var) + epsilon)
            }
        }
    }
    stopifnot(all(FW[FS == 0] == 0))
    ans = list(FW = FW, P = P, losses = na.omit(model$lossess),
               predicates = X$potentials$predicate, call = match.call())
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

check.grad <- function(fun, grad, weight, epsilon = 1e-4) {
    weight = as.double(weight)
    testGrad = grad(weight)
    trueGrad = vector("numeric", length(weight))
    for (i in 1:length(weight)) {
        weight[i] = weight[i] - epsilon
        v1 = fun(weight)
        weight[i] = weight[i] + 2 * epsilon
        v2 = fun(weight)
        trueGrad[i] = (v2 - v1) / epsilon / 2
        c = all.equal(trueGrad[i], testGrad[i], tolerance = epsilon)
        if (!isTRUE(c)) {
            s = c(trueGrad[i], testGrad[i])
            names(s) = c("True", "Test")
            print(s)
            warning(paste0(as.character(c), " at index ", i))
        }
        weight[i] = weight[i] - epsilon
    }
}

logsumexp <- function(x) {
    a = max(x)
    a + log(sum(exp(x - a)))
}

test.me.1 <- function() {
    X = list(c("A"), c("B"), c("A", "B"))
    Y = c(T, F, F)
    fxy.obj <<- fxy(X, Y)
    print(fxy.obj)
    me.obj <<- me(fxy.obj, verbose = FALSE)
    print(me.obj)
}

test.me.2 <- function() {
    spam.sms = na.omit(spam.sms)
    X <<- tk(spam.sms$body)
    Y <<- spam.sms$is.spam
    fxy.obj <<- fxy(X, Y, lower.f = 3L)
    print(fxy.obj)
    me.obj <<- me(fxy.obj, verbose = FALSE)
    print(me.obj)
}
