% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eqMI.semtest.R
\name{eqMI.semtest}
\alias{eqMI.semtest}
\title{Measurement invariance tests using lavaan}
\usage{
eqMI.semtest(..., output = "both", quiet = FALSE)
}
\arguments{
\item{...}{The same arguments as for any lavaan model. See \code{lavaan::cfa} and \code{lavaan::lavOptions} for more information about the arguments.

Users must explicitly specify the name of the input elements for this function to catch. For example, specify 'eqMI.semtest(model = semmodel, data = HolzingerSwineford)' instead just 'eqMI.semtest(semmodel, HolzingerSwineford)'.}

\item{output}{If the function prints out results of covariance structure, mean structure, or both. The value of \code{output} must be \code{mean}, \code{covariance}, or \code{both}. When the tests involve mean structure (\code{output = 'mean' or 'both'}), both the strong and the strict tests of measurement invariance will be conducted.}

\item{quiet}{If quiet=FALSE (default), a summary is printed out containing an overview of the different models that are fitted, together with some model comparison tests. If quiet=TRUE, no summary is printed but results will be stored in the object.}
}
\value{
A list is returned with:
\describe{
\item{\code{LavaanOut}}{A sublist in convention.sem. Contains lavann style output and results for each chi-squared and chi-squared difference.}
\item{\code{Mean.part}}{A sublist in convention.sem. Contains test statistics and fit measures on invariance tests of mean structure.}
\item{\code{Cov.part}}{A sublist in convention.sem. Contains test statistics and fit measures on invariance tests of covariacne structure.}
}
}
\description{
Conventional multiple-group SEM to test measurement invariance. A sequence of chi-squared and chi-squared difference tests will be conducted.
}
\details{
This is a wrapper around the (now deprecated) \code{measurementInvariance} in package \code{semTools}, with the following default options: \code{std.lv = FALSE, fit.measures = "default"}, and \code{method = "satorra.bentler.2001"}. See \code{semTools} for more information. This function is now updated to \code{semTools::measEq.syntax}.
}
\examples{
data(HolzingerSwineford)
semmodel<-'
L1 =~ V1 + V2 + V3
L2 =~ V4 + V5 + V6
L3 =~ V7 + V8
L4 =~ V9 + V10 + V11
'
\donttest{
run.sem <- eqMI.semtest(model = semmodel, data = HolzingerSwineford,
           group = "school", meanstructure = TRUE)
}
}
\references{
Yuan, K. H., & Chan, W. (2016). Measurement invariance via multigroup SEM: Issues and solutions with chi-square-difference tests. Psychological methods, 21(3), 405-426.

Yves Rosseel (2012). lavaan: An R Package for Structural Equation Modeling. Journal of Statistical Software, 48(2), 1-36.

semTools Contributors. (2016). semTools: Useful tools for structural equation modeling. R package version 0.4-14. Retrieved from https://CRAN.R-project.org/package=semTools
}
\seealso{
\code{\link[lavaan]{sem}}, \code{\link[semTools]{measurementInvariance}}
}
\author{
The maintainer, Ge Jiang, adapted the original source code of measurementInvariance() in the \pkg{lavaan} and \pkg{semTools} packages written by Yves Rosseel, Sunthud Pornprasertmanit, and Terrence D. Jorgensen (permission obtained).
}
