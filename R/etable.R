etable <- function(true, guess, ...) {
    stopifnot(is.atomic(true) && is.atomic(guess) && length(true) == length(guess))
    ans = list()
    nx = length(true)
    tb = table(true, guess)
    levels.Y = factor(rownames(tb))
    nc = length(levels.Y)
    ans$summary.table = summary(tb)
    tb = as.matrix(table(true, guess, ...))
    sum.diag = sum(diag(tb))
    tb = rbind(tb, colSums(tb))
    tb = cbind(tb, rowSums(tb))
    tb[nc + 1L, nc + 1L] = sum.diag
    colnames(tb)[nc + 1L] = "Total"
    rownames(tb)[nc + 1L] = "Total"
    ans$accuracy = sum(true == guess) / nx
    ans$table = tb
    type.measure = c("Precision", "Recall", "F-measure", "FPR", "TPR", "AUC")
    measure = matrix(nrow = nc, ncol = length(type.measure))
    colnames(measure) = type.measure
    rownames(measure) = levels.Y
    for (i in 1:nc) {
        target = levels.Y[i]
        tp = sum(true == target & guess == target)
        fp = sum(true != target & guess == target)
        fn = sum(true == target & guess != target)
        tn = sum(true != target & guess != target)
        measure[target, "Precision"] = tp / (tp + fp)
        measure[target, "Recall"] = tp / (tp + fn)
        measure[target, "F-measure"] = 2 / (1 / measure[target, "Precision"] + 1 / measure[target, "Recall"])
    }
    measure = rbind(measure, colMeans(measure))
    rownames(measure)[nc + 1L] = "Average"
    ans$measure = measure
    class(ans) = "etable"
    ans
}


print.etable <- function(obj) {
    print(obj$table)
    print(obj$measure)
    invisible(print(obj$summary.table))
}
