\name{PLEMT}
\alias{PLEMT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Pseudolikelihood based EM test (PLEMT) 
}
\description{
The function conducts the pseudolikelihood based EM test for homogeneity in generalized exponential tilt mixture models
}
\usage{
PLEMT(Y1, Y2, cc = 2, niter = 3, distn = c("beta", "norm", "gamma", "mixture"))
}
%- maybe also 'usage' for other objects documented here.
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
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Chuan Hong, Yong Chen, Yang Ning, Hao Wu 
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
myresult=PLEMT1(dat.case[2,], dat.ctrl[2,], distn="beta")
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

