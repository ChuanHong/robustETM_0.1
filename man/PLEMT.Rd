\name{PLEMT}
\alias{PLEMT}
\title{Pseudolikelihood based EM test (PLEMT)}
\description{
   Conducts the pseudolikelihood based EM test for homogeneity in generalized    exponential tilt mixture models.
}
\usage{
PLEMT(Y1, Y2, cc = 2, niter = 3, distn = c("beta", "norm", "gamma", "mixture"))
}
\arguments{
  \item{Y1}{
the vector of methylation beta values in the case group for one marker
}
  \item{Y2}{
the vector of methylation beta values in the control group for one marker
}
  \item{cc}{
the penalty for the modified pseudolikelihood
}
  \item{niter}{
the number iterations for the EM algorithm
}
  \item{distn}{
name of the distribution for the methylation beta values
}
}

\references{
Submitted to Biometrics:
PLEMT: A Novel Pseudolikelihood Based EM Test for Homogeneity in
Generalized Exponential Tilt Mixture Models
}
\author{
Chuan Hong, Yong Chen, Yang Ning, Hao Wu 
}

\examples{
myresult=PLEMT(dat.case[2,], dat.ctrl[2,], distn="beta")
}

## The function is currently defined as
function (Y1, Y2, cc = 2, niter = 3, distn = c("beta", "norm", 
    "gamma", "mixture")) 
{
    distn <- match.arg(distn)
    n0 <- nrow(Y1)
    result <- rep(0, n0)
    for (i in 1:n0) result[i] <- PLEMT1(Y1[i, ], Y2[i, ], cc, 
        ninter, distn)
    pval <- pnorm(result)
    FDR <- p.adjust(pval, method = "BH")
    res <- data.frame(PLEMT = result, pvalue = , FDR = FDR)
    return(res)
  }
}
\keyword{Conditional likelihood}
\keyword{Penalized likelihood}
\keyword{Semiparametric mixture model}

