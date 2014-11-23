## the main function to perform PLEMT test


## function for performing PLEMT for "one data".
## It takes data from two groups and some configurations, and return test statistics.

PLEMT1 <- function(y1, y2, cc=2, niter=3, distn) {

    logit = function(p) {log(p/(1-p))}
    if(distn=="beta") {
        h1 = function(a) log(a); h2 = function(a) log(1-a)
    } else if (distn=="norm") {
        h1 = function(a) a; h2 = function(a) a^2
    } else if (distn=="mixture") {
        h1 = function(a) a; h2 = function(a) a^2
    } else if (distn=="gamma") {
        h1 = function(a) log(a); h2 = function(a) a
    }

    z = c(y1,y2)
    h1y1 = h1(y1); h1y2 = h1(y2)
    h2y1 = h2(y1); h2y2 = h2(y2)
    n = c(length(y1), length(y2))
    nn = sum(n)

    ## negative log pseudolikelihood
    nploglik <- function(mypar){
        res <- .C("mylikC", as.integer(n[1]), as.integer(n[2]),as.double(h1y1),as.double(h2y1),
                  as.double(h1y2), as.double(h2y2), as.double(mypar),result=double(1))
        return(-res[["result"]]/nn)
    }
    ## negative modified log pseudolikelihood
    nmploglik <- function(t) {
        penalty=cc*log(t[4])
        return(nploglik(t)-penalty)
    }

    beta.init=c(0,0)
    lamda_length=10;alpha_length=10
    lamda_ini=seq(0,1,length=lamda_length)
    alpha_ini=seq(-1,1,length=alpha_length)
    lrt_ini=rep(0,lamda_length)
    TS=hat_lamda=hat_alpha=array(NA, c(lamda_length-1, alpha_length))
    hat_beta=array(NA, c(2, lamda_length-1, alpha_length))

    for (i in 2:lamda_length){
        for (j in 1:alpha_length){
            lamda.temp=lamda_ini[i]
            alpha.temp=alpha_ini[j]
            for (k in 1:niter){
                my.lamda=lamda.temp
                my.alpha=alpha.temp
                nmploglik_EM_beta=function(t, alpha, lamda){
                    nmploglik(c(alpha,t,lamda))
                }
                op <- optim(beta.init, nmploglik_EM_beta, alpha=my.alpha, lamda=my.lamda,
                            hessian=TRUE, method = "Nelder-Mead", control = list(fnscale=1,maxit=1000))
                my.beta=op$par
                upper=my.lamda*exp(my.alpha+my.beta[1]*h1y2+my.beta[2]*h2y2)
                lower=1-my.lamda+my.lamda*exp(my.alpha+my.beta[1]*h1y2+my.beta[2]*h2y2)
                w=upper/lower
                if(any(is.na(w))==FALSE){
                    my.lamda=(sum(w)+1)/(length(y2)+1)
                    op <- optim(beta.init, nmploglik_EM_beta, alpha=my.alpha,lamda=my.lamda,
                                hessian=TRUE, method = "Nelder-Mead", control = list(fnscale=1,maxit=1000))
                    my.beta=op$par
                    alpha.temp=my.alpha
                    lamda.temp=my.lamda}

                if(any(is.na(w))==TRUE){
                    my.beta=my.beta
                    alpha.temp=my.alpha
                    lamda.temp=my.lamda}
            }

            hat_lamda[i-1, j]=my.lamda
            hat_alpha[i-1, j]=my.alpha
            hat_beta[,i-1, j]=my.beta
            TS[i-1, j]=2*(-nmploglik(c(my.alpha, my.beta, my.lamda))+nmploglik(c(0,0,0,1)))
        }
    }

    mplrt_EM.TS=2*max(TS)
    return(mplrt_EM.TS)
}


## A wrapper function, takes two matrices with same number of rows.
## Each row is for a "gene" or CpG site. This returns other things like p-values and FDR.
PLEMT <- function(Y1, Y2, cc=2, niter=3, distn=c("beta","norm", "gamma", "mixture")) {
    distn <- match.arg(distn)

    ## loop over CG sites
    n0 <- nrow(Y1)
    result <- rep(0, n0)
    for(i in 1:n0)
        result[i] <- PLEMT1(Y1[i,], Y2[i,], cc, niter, distn)

    ## make other fields, like p-values, FDR, etc.
    ## Modify the following two lines for computing p-values and FDR, if needed.
    pval <- pchisq(result, df=2, lower.tail=FALSE)
    FDR <- p.adjust(pval, method="BH")
    res <- data.frame(PLEMT=result, pvalue=pval, FDR=FDR)
    return(res)
}

