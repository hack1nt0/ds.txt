me <- function(X, lambda = 0, alpha = 0, method = c("Adam", "RAdam", "L-BFGS-B"),
               maxit = 100L, tolerance = 1e-9, eta = 0.001, beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8, AMSGrad = FALSE, decay = FALSE, verbose = FALSE, keep = FALSE,
               FW = matrix(rnorm(nf), nrow = np, ncol = nc)) {
    stopifnot(is.fxy(X))
    ans = list()
    ans$fxy.obj = X
    ans$levels.Y = X$levels.Y
    env <- new.env(parent = emptyenv())
    env$losses = rep(NA, maxit)
    env$it = 0L
    Y = X$Y
    FX = X$FX
    FS = X$FS
    nx = nrow(X$FX)
    nc = ncol(X$FS)
    np = nrow(X$FS)
    nf = np * nc
    FW[FS == 0] = 0
    method = match.arg(method)
    if (method == "Adam") {
        tFX = t(FX)
        res = .Call("C_me", PACKAGE = "ds.txt", nx, np, nc, tFX@p, tFX@i, FX@p, FX@i, Y, lambda, alpha, FW, FS, maxit,
                    tolerance, eta, beta1, beta2, epsilon, AMSGrad, decay, verbose, keep)
        ans$W = FW
        if (keep) {
            ans$P = res$P
            ans$mean.gr = res$mean.gr
            ans$var.gr = res$var.gr
            ans$max.var.gr = res$max.var.gr
            ans$losses = res$losses
        } else {
            ans$losses = res
        }
    }
    else {
        ixs = 1:nx; ips = 1:np
        adjDfm = adjacent.list(FX, byrow = T)
        adjFdm = adjacent.list(FX, byrow = F)
        FG = matrix(0, nrow = np, ncol = nc)
        P = matrix(0, nrow = nx, ncol = nc)
        FC = matrix(0, nrow = np, ncol = nc)
        for (ip in 1:np) {
            ys = Y[adjFdm[[ip]]]
            for (y in ys)
                FC[ip, y] = FC[ip, y] + 1
        }
        loss <- function(weight) {
            FW = matrix(weight, nrow = np, ncol = nc)
            nx = length(adjDfm)
            lhs = 0
            rhs = 0
            for (ix in 1:nx) {
                ifs = adjDfm[[ix]]
                lse = logsumexp(Matrix::colSums(FS[ifs,, drop = F] * FW[ifs,, drop = F])) #todo
                lhs = lhs + lse
                rhs = rhs + sum(FS[ifs, Y[ix]] * FW[ifs, Y[ix]])
            }
            curLoss = (lhs - rhs) / nx + lambda * ((1 - alpha) * sum(FW ^ 2) + alpha * sum(abs(FW))) #todo
            # curLoss = (lhs - rhs) / nx
            if (verbose)
                print(paste0("it ", env$it, ", loss ", curLoss))
            env$losses[env$it] = curLoss
            env$it = env$it + 1L
            curLoss
        }
        grad <- function(FW) {
            FW = matrix(FW, nrow = np, ncol = nc)
            for (ix in ixs) {
                ifs = adjDfm[[ix]]
                P[ix,] = Matrix::colSums(FS[ifs,, drop = F] * FW[ifs,, drop = F]) #todo
                P[ix,] = exp(P[ix,] - logsumexp(P[ix,]))
            }
            # P[,] = exp(P)
            # P[,] = P / Matrix::rowSums(P)
            # P[,] = P / rowSums(P) #todo necessary?
            for (ip in ips) {
                ids = adjFdm[[ip]]
                FG[ip,] = Matrix::colSums(P[ids,, drop = F])
            }
            FG[,] = (FG - FC) * FS / nx + lambda * ((1 - alpha) * 2 * FW + alpha * sign(FW))
            as.double(FG) # todo
        }
        # check.grad(loss, grad, FW, epsilon = 1e-9)
        if (method == "RAdam") {
            gr = matrix(0, np, nc)
            gr1 = matrix(0, np, nc)
            gr2 = matrix(0, np, nc)
            max.gr2 = matrix(0, np, nc)
            powBeta1 = beta1
            powBeta2 = beta2
            for (it in 1:maxit) {
                if (verbose)
                    loss(FW)
                gr[,] = grad(FW)
                gr1[,] = beta1 * gr1 + (1 - beta1) * gr
                gr2[,] = if (decay) (powBeta2 * gr2 + (1 - powBeta2) * gr * gr)  else (beta2 * gr2 + (1 - beta2) * gr * gr)
                if (AMSGrad) {
                    max.gr2[,] = pmax(max.gr2, gr2)
                    FW[,] = FW - eta * gr1 / (sqrt(max.gr2) + epsilon)
                }
                else
                    FW[,] = FW - eta * (gr1 / (1 - powBeta1)) / (sqrt(gr2 / (1 - powBeta2)) + epsilon)
                powBeta1 = powBeta1 * beta1
                powBeta2 = powBeta2 * beta2
            }
            ans$W = FW
        }
        else if (method == "L-BFGS-B") {
            opt.obj = optim(as.double(FW), fn = loss, gr = grad, method = "L-BFGS", control = list(maxit = maxit))
            ans$W = matrix(opt.obj$par, nrow = np, ncol = nc)
        }
        ans$losses = na.omit(env$losses)
    }
    stopifnot(all(FW[FS == 0] == 0))
    ans$predicates = X$potentials$predicate
    ans$names = X$potentials$name
    ans$call = match.call()
    class(ans) = "me"
    ans
}

cv.me <- function(X, lambdas = NULL, alphas = c(0, 0.5, 1),
                  type.measure = c("class", "auc", "log.likelihood"),
                  nfold = 10L, foldid = NULL, parallel = FALSE,
                  maxit = 100L, tolerance = 1e-9, eta = 0.001, beta1 = 0.9, beta2 = 0.999, epsilon = 1e-8, AMSGrad = FALSE, verbose = FALSE,
                  FW = matrix(rnorm(nf), nrow = np, ncol = nc)) {
    type.measure = match.arg(type.measure)
    # type.predict = switch(type.measure, class = "class", auc = "response", log.likelihood = "response")
    lambdas = unique(sort(lambdas, decreasing = T))
    nlambda = length(lambdas)
    alphas = unique(alphas)
    nalpha = length(alphas)
    FX = X$FX
    tFX = t(FX)
    FS = X$FS
    Y = X$Y
    nx = nrow(FX)
    np = nrow(FS)
    nc = ncol(FS)
    nf = np * nc
    nfold = as.integer(min(nx, nfold))
    batches = NULL
    if (is.null(foldid)) {
        batches = make.batches(nx, nbatch = nfold)
    } else {
        foldid = as.integer(foldid)
        min.foldid = min(foldid)
        max.foldid = max(foldid)
        stopifnot(1 <= min.foldid && max.foldid <= nfold)
        batches = adjacent.list(foldid, byrow = F, ngroup = nfold)
        stopifnot(all(sapply(batches, length) > 0L))
    }
    grids = expand.grid(lambda = lambdas, alpha = alphas)
    losses.batches = vector(mode = "list", length = nfold)
    for (foldi in 1:nfold) {
        batch = batches[[foldi]]
        train.X = FX[-batch,]
        t.train.X = tFX[-batch,]
        train.Y = Y[-batch,]
        t.test.X = tFX[batch,]
        test.Y = Y[batch,]
        res = .Call("C_cv_me", PACKAGE = "ds.txt", t.train.X@p, t.train.X@i, train.X@p, train.X@i, train.Y, t.test.X@p, t.test.X@i, test.Y, lambdas, alphas, type.measure, FW, FS, maxit, tolerance, eta, beta1, beta2, epsilon, AMSGrad, verbose)
        losses.batches[[foldi]] = res
    }
    losses = do.call(`+`, losses.batches) / nfold
    ans = cbind(grids, losses)
    class(ans) = "cv.me"
    ans
}

make.batches <- function(x, nbatch = 10L, random = FALSE) {
    stopifnot(is.integer(x) && all(x > 0) && is.integer(nbatch) && nbatch > 0 && is.logical(random))
    ans = vector("list", nbatch)
    if (length(x) == 1L)
        x = 1:x
    n = length(x)
    nbatch = min(nbatch, n)
    if (random) x = sample(x, size = n)
    size.batch = as.integer((n + nbatch - 1) / nbatch)
    for (i in 1:nbatch)
        ans[[i]] = x[((i - 1L) * size.batch + 1L) : min(i * size.batch, n)]
    ans
}

predict.me <- function(me.obj, X, type = c("class", "response", "link", "coefficients", "nonzero"),
                       target.class = me.obj$levels.Y[[1L]], threshold = 0.5, ...) {
    ans = list()
    type = match.arg(type)
    nx = NULL
    levels.Y = me.obj$fxy.obj$levels.Y
    ans$levels.Y = levels.Y
    adjDfm = NULL
    if (is.fxy(X)) {
        adjDfm = adjacent.list(X$FX, byrow = T)
        nx = nrow(X$FX)
    }
    else {
        FX = fx(X, predicates = me.obj$predicates, names = me.obj$names)$fx
        adjDfm = adjacent.list(FX, byrow = T, dict = colnames(me.obj$fxy.obj$FX)) #todo
        nx = length(X)
    }
    nc = length(levels.Y)
    W = me.obj$W
    S = me.obj$fxy.obj$FS
    np = nrow(W)
    if (type == "class" || type == "response" || type == "link") {
        P = matrix(0, nx, nc)
        colnames(P) = levels.Y
        for (ix in 1:nx) {
            ifs = adjDfm[[ix]]
            P[ix,] = Matrix::colSums(W[ifs,, drop = F] * S[ifs,, drop = F])
            if (type != "link")
                P[ix,] = exp(P[ix,] - logsumexp(P[ix,]))
        }
        if (type == "class") {
            guess.Y = NULL
            if (nc == 2L) {
                complement.class = NULL
                if (levels.Y[1L] == target.class)
                    complement.class = levels.Y[2L]
                else
                    complement.class = levels.Y[1L]
                print(complement.class)
                guess.Y = rep(complement.class, nx)
                guess.Y[P[, target.class] >= threshold] = target.class
            } else {
                guess.Y = levels.Y[apply(P, 1L, which.max)]
            }
            ans$result = guess.Y
        } else
            ans$result = P
    }
    else {
        # Type(s) for debug or interactive coding.
        stopifnot(!is.fxy(X))
        result = vector(mode = "list", length = nx)
        name.predicates = colnames(me.obj$fxy.obj$FX)
        for (ix in 1:nx) {
            ifs = adjDfm[[ix]]
            evidence = W[ifs,, drop = F] * S[ifs,, drop = F]
            if (type == "nonzero") {
                rownames(evidence) = name.predicates[ifs]
                colnames(evidence) = levels.Y
                result[[ix]] = evidence
            } else {
                xi = fx(X[[ix]], predicates = me.obj$predicates, names = me.obj$names)$fx
                name.predicates.2 = colnames(xi)
                evidence2 = Matrix(0, nrow = length(name.predicates.2), ncol = nc)
                rownames(evidence2) = name.predicates.2
                colnames(evidence2) = levels.Y
                evidence2[name.predicates[ifs],] = evidence
                result[[ix]] = evidence2
            }
        }
        ans$result = result
        ans$X = X
    }
    ans$type = type
    class(ans) = "predict.me"
    ans
}

print.predict.me <- function(obj, order.by = NULL) {
    type = obj$type
    if (type == "nonzero" || type == "coefficients") {
        if (!is.null(order.by))
            order.by = match.arg(order.by, obj$levels.Y)
        evidences = obj$result
        nx = length(evidences)
        for (i in 1:nx) {
            ei = evidences[[i]]
            M = rbind(ei, colSums(ei))
            rownames(M)[nrow(M)] = "Total:"
            M = as(M, "dgCMatrix")
            if (!is.null(order.by)) {
                ord = order(M[, order.by])
                M[,] = M[ord,]
            }
            if (i > 1L)
                cat("---------------------------------------\n")
            cat("Text: ");
            xi = obj$X[[i]]
            if (length(xi) == 1L)
                cat(xi, "\n")
            else
                print(xi)
            cat("Evidence: "); invisible(print(M));
        }
    } else {
        warning("Not implemented yet.")
    }
}

plot.me <- function(obj, ...) {
    plot(obj$train.losses, xlab = "Iteration", ylab = "Negative Log Likelihood", t = 'o', ...)
    points(obj$valid.losses, col = "red")
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
    set.seed(1)
    fxy.obj <- fxy(c("A", "B", "AB"), c(T, F, T))
    print(fxy.obj)
    me.obj <- me(fxy.obj, method = "Adam", maxit = 10L, verbose = T)
    print(me.obj)
    # pred.me <- predict(me.obj, c("C B", "A"), type = "coef")
    # print(pred.me, order.by = "T")
}

test.me.3 <- function() {
    fxy.obj <- fxy(c("A", "B", "AB"), c(T, F, T))
    print(fxy.obj)
    cv.me.obj <- cv.me(fxy.obj, maxit = 1L, verbose = T)
    cv.me.obj
}

test.me <- function(...) {
    set.seed(1)
    spam.sms = na.omit(spam.sms)
    nx = nrow(spam.sms)
    train.idx = sample.int(nx, size = as.integer(3 / 4 * nx), replace = F)
    train.X = spam.sms$body[train.idx]
    train.Y = spam.sms$is.spam[train.idx]
    test.X = spam.sms$body[-train.idx]
    test.Y = spam.sms$is.spam[-train.idx]
    fxy.obj <<- fxy(train.X, train.Y, lower.f = 3L)
    me.obj <<- me(fxy.obj, ..., verbose = T)
    pred.me <<- predict(me.obj, test.X, type = "class")
    print(etable(test.Y, pred.me$result))
}
