\name{robustETM-package}
\alias{robustETM-package}
\alias{robustETM}
\docType{package}
\title{
Robust tests for exponential tilt models (robustETM)
}
\description{
Conduct robust tests for homogeientiy in exponential tilt models
}
\details{
\tabular{ll}{
Package: \tab robustETM\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-11-22\cr
License: GPL (version 2 or later)
}
}
\author{
Chuan Hong, Yong Chen, Yang Ning, Hao Wu
Maintainer: Chuan Hong <hong.chuan.hannah@gmail.com>
}
\references{
Submitted to Biometrics

PLEMT: A Novel Pseudolikelihood Based EM Test for Homogeneity in
Generalized Exponential Tilt Mixture Models
}
\examples{
data(pseudo_dat)
myresult=PLEMT(dat.case[c(1:2),], dat.ctrl[c(1:2),], cc=2, niter=3, distn="beta")
}
\keyword{robustETM}
\keyword{PLEMT}
\keyword{Conditional likelihood}
\keyword{Penalized likelihood}
\keyword{Semiparametric mixture model}

