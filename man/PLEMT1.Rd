\name{PLEMT1}
\alias{PLEMT1}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
PLEMT1
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
PLEMT1(y1, y2, cc = 2, niter = 3, distn)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y1}{
%%     ~~Describe \code{y1} here~~
}
  \item{y2}{
%%     ~~Describe \code{y2} here~~
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
function (y1, y2, cc = 2, niter = 3, distn) 
{
    logit = function(p) {
        log(p/(1 - p))
    }
    if (distn == "beta") {
        h1 = function(a) log(a)
        h2 = function(a) log(1 - a)
    }
    else if (distn == "norm") {
        h1 = function(a) a
        h2 = function(a) a^2
    }
    else if (distn == "mixture") {
        h1 = function(a) a
        h2 = function(a) a^2
    }
    else if (distn == "gamma") {
        h1 = function(a) log(a)
        h2 = function(a) a
    }
    z = c(y1, y2)
    h1y1 = h1(y1)
    h1y2 = h1(y2)
    h2y1 = h2(y1)
    h2y2 = h2(y2)
    n = c(length(y1), length(y2))
    nn = sum(n)
    nploglik <- function(mypar) {
        res <- .C("mylikC", as.integer(n[1]), as.integer(n[2]), 
            as.double(h1y1), as.double(h2y1), as.double(h1y2), 
            as.double(h2y2), as.double(mypar), result = double(1))
        return(-res[["result"]]/nn)
    }
    nmploglik <- function(t) {
        penalty = cc * log(t[4])
        return(nploglik(t) - penalty)
    }
    beta.init = c(0, 0)
    lamda_length = 10
    alpha_length = 10
    lamda_ini = seq(0, 1, length = lamda_length)
    alpha_ini = seq(-1, 1, length = alpha_length)
    lrt_ini = rep(0, lamda_length)
    TS = hat_lamda = hat_alpha = array(NA, c(lamda_length - 1, 
        alpha_length))
    hat_beta = array(NA, c(2, lamda_length - 1, alpha_length))
    for (i in 2:lamda_length) {
        for (j in 1:alpha_length) {
            lamda.temp = lamda_ini[i]
            alpha.temp = alpha_ini[j]
            for (k in 1:niter) {
                my.lamda = lamda.temp
                my.alpha = alpha.temp
                nmploglik_EM_beta = function(t, alpha, lamda) {
                  nmploglik(c(alpha, t, lamda))
                }
                op <- optim(beta.init, nmploglik_EM_beta, alpha = my.alpha, 
                  lamda = my.lamda, hessian = TRUE, method = "Nelder-Mead", 
                  control = list(fnscale = 1, maxit = 1000))
                my.beta = op$par
                upper = my.lamda * exp(my.alpha + my.beta[1] * 
                  h1y2 + my.beta[2] * h2y2)
                lower = 1 - my.lamda + my.lamda * exp(my.alpha + 
                  my.beta[1] * h1y2 + my.beta[2] * h2y2)
                w = upper/lower
                if (any(is.na(w)) == FALSE) {
                  my.lamda = (sum(w) + 1)/(length(y2) + 1)
                  op <- optim(beta.init, nmploglik_EM_beta, alpha = my.alpha, 
                    lamda = my.lamda, hessian = TRUE, method = "Nelder-Mead", 
                    control = list(fnscale = 1, maxit = 1000))
                  my.beta = op$par
                  alpha.temp = my.alpha
                  lamda.temp = my.lamda
                }
                if (any(is.na(w)) == TRUE) {
                  my.beta = my.beta
                  alpha.temp = my.alpha
                  lamda.temp = my.lamda
                }
            }
            hat_lamda[i - 1, j] = my.lamda
            hat_alpha[i - 1, j] = my.alpha
            hat_beta[, i - 1, j] = my.beta
            TS[i - 1, j] = 2 * (-nmploglik(c(my.alpha, my.beta, 
                my.lamda)) + nmploglik(c(0, 0, 0, 1)))
        }
    }
    mplrt_EM.TS = 2 * max(TS)
    return(mplrt_EM.TS)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
