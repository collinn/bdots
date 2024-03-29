% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bdotsFit.R
\name{bdotsFit}
\alias{bdotsFit}
\title{Fit nlme curves to grouped observations}
\usage{
bdotsFit(
  data,
  subject,
  time,
  y,
  group,
  curveType = doubleGauss(concave = TRUE),
  cor = TRUE,
  numRefits = 0,
  cores = 0,
  verbose = FALSE,
  ...
)
}
\arguments{
\item{data}{Dataset used}

\item{subject}{Column name of dataset containing subject identifiers}

\item{time}{Column name containing time variable}

\item{y}{Column name containing outcome of interest}

\item{group}{Character vector containing column names of groups. Can be
greater than one}

\item{curveType}{See details/vignette}

\item{cor}{Boolean. Autocorrelation?}

\item{numRefits}{Integer indicating number of attempts to fit an observation
if the first attempt fails}

\item{cores}{number of cores. Default is \code{0}, indicating half cores available}

\item{verbose}{currently not used}

\item{...}{Secret}
}
\value{
Object of class 'bdotsObj', inherits from data.table
}
\description{
Creates observation level curves to use in bdotsBoot
}
\details{
This is step one of the three step bdots process. Things should be
more or less straight forward. The only tricky part involves curveType. For now
know that one can use doubleGauss(concave = TRUE/FALSE) or logistic(). Should
be passed in as a call. See the vignette on customizing this
}
\examples{
\dontrun{
res <- bdotsFit(data = cohort_unrelated,
                subject = "Subject",
                time = "Time",
                y = "Fixations",
                group = c("Group", "LookType"),
                curveType = doubleGauss(concave = TRUE),
                cor = TRUE,
                numRefits = 2,
                cores = 0,
                verbose = FALSE)
}

}
