library(data.table)
library(bdots)

## fit curves with Group and Looktype
fit <- bdotsFit(data = cohort_unrelated,
                subject = "Subject",
                time = "Time",
                y = "Fixations",
                group = c("Group", "LookType"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 2,
                cores = 0,
                verbose = FALSE)


## using data.table::fread
mat <- fread("~/tmp/AUDIO_7_Cohort_out.txt")

## Only need subject and group columns + parameters
mat <- mat[, 1:8]

# to match LookType to values
mat$LookType <- ifelse(mat$Group == "C", "Cohort", "Unrelated_Cohort")

# Group (values 50, 65) not provided, so pretending they were 50
mat$Group <- 50


bdo <- fit
x <- mat
mergeMatlabFits <- function(x, bdo) {

  ## merge columns

}
