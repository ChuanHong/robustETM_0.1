\name{PLEMT}
\alias{PLEMT}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PLEMT(Y1, Y2, cc = 2, niter = 3, distn = c("beta", "norm", "gamma", "mixture"))
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Y1}{
%%     ~~Describe \code{Y1} here~~
}
  \item{Y2}{
%%     ~~Describe \code{Y2} here~~
}
  \item{cc}{
%%     ~~Describe \code{cc} here~~
}
  \item{niter}{
%%     ~~Describe \code{niter} here~~
}
  \item{distn}{
%%     ~~Describe \code{distn} here~~
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
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (Y1, Y2, cc = 2, niter = 3, distn = c("beta", "norm", 
    "gamma", "mixture")) 
{
    distn <- match.arg(distn)
    n0 <- nrow(Y1)
    result <- rep(0, n0)
    for (i in 1:n0) result[i] <- PLEMT1(Y1[i, ], Y2[i, ], cc, 
        ninter, distn)
    pval <- pchisq(result, df = 2, lower.tail = FALSE)
    FDR <- p.adjust(pval, method = "BH")
    res <- data.frame(PLEMT = result, pvalue = pval, FDR = FDR)
    return(res)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
