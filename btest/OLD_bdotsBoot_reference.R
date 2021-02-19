bdotsBoot <- function(ff, bdObj, N.iter = 1000, alpha = 0.05, p.adj = "oleson", cores = 0, ...) {

    if (cores < 1) cores <- detectCores()/2

    ## Could maybe list what was removed
    if(any(bdObj$fitCode == 6)) {
        warning("Some observations had NULL gnls fits. These will be removed")
        bdObj <- bdObj[fitCode != 6, ]
    }

    ## pretend ff is parsed. This gives us diffGroup, compareGroup
    #diffGroup <- "LookType" # or NULL
    #compareGroup <- "Group"

    #, ah, NOPE!
    ## Split by subject and diffgroup. Each resulting
    ## diffGroup can be null, and this splits correctly
    #test <- split(bdObj, by = c("Subject", NULL), drop = TRUE)
    ## I think we need to go a level up. First  split by compare group, then split
    # each of those by subject
    ## Actually, that affirms we need a name change.  For Bob's data, we have Group (LI/TD)
    # But first, if diffs = TRUE, we split by diffGroup, so that we have a list
    # of 2, and on that list of two, we lapply split by subject
    diffList <- split(bdObj, by = c("Subject", diffGroup), drop = TRUE)
    diffList <- lapply(diffList, function(x) {class(x) <- class(bdObj); x})

    ##  First split by diffGroup (could be NULL). Then get corMat if necessary
    ## mm....should probably determine if this is paired now
    # which, actually would be true if each item if diffList is nrow2
    ## Determine if paired
    ## Should not need to error check this. If every split for diffList is length
    # 2, then paired. cor would only fail, then, if number of pars were different
    is.paired <- !any(vapply(diffList, nrow, numeric(1)) != 2) # this is wrong
    if (is.paired) {
        coefs <- lapply(split(bdObj, by = compareGroup, drop = TRUE), function(x) {
            class(x) <- class(bdObj)
            coef.bdots(x)
        })
        corMat <- do.call(cor, setNames(coefs, c("x", "y")))
    } else {
        corMat <- NULL
    }


    ## Draw the bootstraps
    ## res is a length diffList list of matrices for each subject
    ## NOTE! THIS DOES NOT DIFFERENTIATE COMPARE GROUP
    # which will be N.iter x (numpars (x2 if paired))
    res <- lapply(diffList, bdotsBooter, N.iter, corMat)


    ## Am I taking the same "mean" as Jake intended in his paper?
    # am I applying this the correct way?
    ## This is wrong, it's taking mean across Group = LI and Group = TD in Bob's data
    # because this is N.iter x numPars(x2 if paired) of the mean for each iteration
    meanMat <- Reduce(`+`, res)/length(res) # maybe no, this is mean across subjects for each iterations. Actually, maybe yes

    ## Get formula and turn into function
    ## For example, a user might write y ~ f(time1, expr) + g(time2, expr) if they wanted
    # to try to fit two curves at once (or by Subj == w/e). That should be  made a priority
    # i.e., my data.frame shouldn't be named time everywhere (or y, for that matter)
    ff <- attr(bdObj, "formula")
    f_bod <- deparse1(ff[[3]])
    f_args <- paste0(colnames(coef.bdots(bdObj)), collapse = ", ") # + colnames(dat)
    eval(parse(text = paste('curveFun <- function(', f_args, ', time', ') { return(' , f_bod , ')}', sep='')))

    ## Not the correct way to do this
    if (is.paired) {

        time <- attr(bdObj, "time")

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
        res1 <- lapply(mmList, function(x){force(x); do.call(curveFun, x)}) # this is 1000 curve fits

        mm <- meanMat[, (ncol(meanMat)/2 + 1):ncol(meanMat)]
        mmList <- split(mm, row(mm))
        mmList <- lapply(mmList, function(x) {
            x <- as.list(x)
            x$time <- time
            setNames(x, c(parNames, "time"))
        })
        res2 <- lapply(mmList, function(x){force(x); do.call(curveFun, x)})

        ### This is, in general, the idea of differencing curves
        ### so we could wrap this up into a single idea as well
        ## This is paired, so next steps
        # i. take mean of each element in the list
        # ii. subtract them
        #   a. these had to be matrices, not vectors first. False. Use `+` not sum (good to know)
        res11 <- Reduce(`+`, res1) / length(res1)
        res22 <- Reduce(`+`, res2) / length(res2)



        ### The below occurs if I leave Subject == 405 in Bob's data
        # However, this method in general is obsolete, so have to logic back
        # around to reproduce it
        ## (names of things may change, I believe here mmList is the split
        ## for the 4 parameters of the  first group)
        ## GOT IT.  SOME OF THE SLOPE VALUES ARE NEGATIVE (FOR BOBS DATA)
        ## res[[7]] has seriously fucked slopes
        # This is it. REMOVING THIS SUBJECT FIXED IT!!! WOOHOO!
        # need to figure out why it was so bad
        # > names(res)[7]
        # [1] "405.LI"
        # rr <- vector("list", length = length(mmList))
        # for(i in  seq_along(mmList)) {
        #   MM <- mmList[[i]]
        #   mini <- MM$mini
        #   peak <- MM$peak
        #   slope <- MM$slope
        #   cross <- MM$cross
        #   TIME <- MM$time
        #   rr[[i]] <- curveFun(mini, peak, slope, cross, TIME)
        # }
        # rr1 <- Reduce(`+`, rr) / length(rr)


        ## These didn't work, but illustrating example of why not, given above
        # mmListprep <- lapply(mmList, function(x) {x$time <- unlist(x$time); x})
        # res <- lapply(mmListprep, function(x) do.call(curveFun, as.list(x))) # <- issue here with time being a list

        ## Ok, let's do it again for 2, but obviously this is not how it will be
        ## its still embarassing to get the idea out. But why? I think this is ok as a
        # temporary solution. And it's a MUCH better way of keeping notes of ideas, really
        # This kind of part of my writing step



    }

}
