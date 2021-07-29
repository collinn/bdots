library(bdots)

groups <- c("Group", "LookType")
res <- bdotsFit(data = cohort_unrelated,
                subject = "Subject",
                time = "Time",
                y = "Fixations",
                group = groups,
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 2,
                cores = 2,
                verbose = FALSE)

### Verify components of data.table

## Verify class correct, and dt info
expect_equal(class(res), c("bdotsObj", "data.table", "data.frame"))

## Correct column
# Subject, y, fit, ar1, r2, fitcode
expect_equal(ncol(res), 5 + length(groups))

## Verify formula is fine
ff <- attr(res, "formula")

expect_equal(class(ff), "call")
expect_equal(names(nn <- attributes(ff)), "parnames")

# par names present in formula?
ss <- Reduce(`+`, lapply(unlist(nn), grep, x = deparse1(ff)))
expect_equal(ss, length(unlist(nn)))

## Grouping var match dataset?
grp_vals <- do.call(function(...) paste(..., sep = "."), unique(cohort_unrelated[, ..groups]))
grp_attr <- attr(res, "groups")

expect_equal(grp_attr[['groups']], groups)
expect_equal(grp_attr[['vals']], grp_vals)
