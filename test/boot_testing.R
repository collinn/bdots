
bdObj <- res.b; diffGroup = "Group"; compareGroup <- "TrialType"

## To get list of cor mats

# each of those by subject and diffGroup
# if diff == TRUE then do this
cmpList <- split(bdObj, by = diffGroup, drop = TRUE) # list of the 2 groups
## Things should be paired if the subject IDs are the same in all components of cmpList
# in other words, this represents the "outer" split
corMat <- lapply(cmpList, function(x) {
    tt <- split(x, by = compareGroup)
    cm <- lapply(tt, coef.bdots)
    do.call(cor, setNames(cm, c("x", "y")))
})

cmpList <- lapply(cmpList, function(x) split(x, by = "Subject")) # nested list of what we had done previously
# for example, cmpList[[1]][[1]] is sub = 148, group = LI, TrialType = M/W

## if paired
time <- attr(bdObj, "time")

## As is, this puts us at length 2 list, with
# each element containing the N.iter parameter samples for each subject
# This should apply for diffs = TRUE. We take curve1 - curve2 for each cmpList group
# If diffs = FALSE, then cmpList would only be length one (should change name, then)
hmm <- lapply(names(cmpList), function(nn) {
    x <- cmpList[[nn]]
    is.paired <- !any(vapply(x, nrow, numeric(1)) != 2)
    # length - number of sujbects in group x
    # element is N.iter x pars
    res <- lapply(x, bdotsBooter, N.iter, corMat[[nn]])
    meanMat <- Reduce(`+`, res)/length(res)


    ## Try this (we already know it's not ideal)
    ## Faster (by FAR) not in parallel
    mm <- meanMat[, 1:(ncol(meanMat)/2)]
    mmList <- split(mm, row(mm)) # should I check on off chance not unique?
    if(length(mmList) != nrow(mm)) stop("crazy anomalie of probability 0. Run again w/ different seed")

    parNames <- colnames(mm)
    mmList <- lapply(mmList, function(x) {
        x <- as.list(x)
        x$time <- time
        setNames(x, c(parNames, "time"))
    })
    ## I should unlist this to get numeric matrix, and then rebuild it
    res1 <- lapply(mmList, function(x){force(x); do.call(curveFun, x)}) # this is 1000 curve fits
    res1 <- matrix(unlist(res1, use.names = FALSE), nrow = length(res1), byrow = TRUE)
    curve1 <- colMeans(res1) # re bobs data, this is curve for TrialType = M
    sd1 <- apply(res1, 2, sd)


    mm <- meanMat[, (ncol(meanMat)/2 + 1):ncol(meanMat)]
    mmList <- split(mm, row(mm))
    mmList <- lapply(mmList, function(x) {
        x <- as.list(x)
        x$time <- time
        setNames(x, c(parNames, "time"))
    })
    res2 <- lapply(mmList, function(x){force(x); do.call(curveFun, x)})
    res2 <- matrix(unlist(res2, use.names = FALSE), nrow = length(res2), byrow = TRUE)
    curve2 <- colMeans(res2) # re bobs data, this is curve for TrialType = W
    sd2 <- apply(res2, 2, sd)
})

