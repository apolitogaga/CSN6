library(igraph)
library(tolerance)
library(xtable)
library(VGAM)
library(stats4)



#-LogLikelyhood functions define
minus_LL_geometric <- function(q) {
    N <- length(x)
    M <- sum(x)
    -(M-N)*log(1-q)-N*log(q)
}

minus_LL_poisson <- function(lambda) {
    N <- length(x)
    M <- sum(x)
    #C <- sum(log(factorial(x)))
    #SUM[1..N](SUM[2..ki] log j)
    C <- sum(sapply(x,function(y) sum(log(2:y))))
    -M*log(lambda)+N*(lambda+log(1-exp(-lambda)))+C
}

minus_LL_zeta <- function(gamma) {
    length(x) * log(zeta(gamma)) + gamma * sum(log(x))
}

minus_LL_zeta_gamma2 <- function() {
    3 * sum(log(x)) + length(x) * log(zeta(3, deriv = 0))
}


minus_LL_zeta_trunc <- function(gamma, Kmax) {
    H <- 0
    for(i in 1:Kmax) {
        H <- H + i^(-gamma)
    }
    gamma * sum(log(x)) + length(x) * log(H)
}

minus_log_like_zetaright <- function(gamma, Kmax) {
    gamma * M_dash + N * log(sum(seq(1,Kmax)^-gamma))
}


#This two functions seem to not work properly, so we are using the functions given by Tolerance package
zipf_LL <- function(s) sum(x*(s*log(1:length(x))+log(sum(1/(1:length(x))^s))))
zipfMan_LL <-  function(s, q) sum(x * (s * log(1:length(x) + q) + log(sum(1/(1:length(x) + q)^s))))



#AIC
get_AIC <- function(m2logL,K,N) {
    m2logL + 2*K*N/(N-K-1) # AIC with a correction for sample size
}

#MLE Estimation
findBestModel <- function(gamma_t0) {
    N <- length(x)
    mu <- sum(x)/N

    #Calculate parameters for each distribution
    mle_zeta_trunc <- tryCatch(
                        mle(minus_LL_zeta_trunc,
                          start = list(Kmax=max(x),gamma=3),
                          method = "L-BFGS-B",
                          lower = c(N, 1),
                          upper = c(Inf,Inf)),
                        error=function(e) NULL)


    mle_zeta <- mle(minus_LL_zeta,
                    start = list(gamma=3),
                    method = "L-BFGS-B",
                    lower = c(1.000001))

    mle_poisson <-  mle(minus_LL_poisson,
                        start = list(lambda=mu),
                        method = "L-BFGS-B",
                        lower = c(0.000001))

    mle_geometric <- mle(minus_LL_geometric,
                         start = list(q=(1/mu)),
                         method = "L-BFGS-B",
                         lower = c(0.0000001),
                         upper = 0.999999)

    #mle_zipf <- mle(zipf_LL,
    #           start=list(s=1),
    #           method="L-BFGS-B",
    #           lower=c(1.0000001))
    # mle_zipf <- zm.ll(x, dist="Zipf", s=1)

    #mle_zipfman <- mle(zipfMan_LL,
    #                start=list(s=1,q=1),
    #                method="L-BFGS-B",
    #                lower=c(1,0.0001))

    # mle_zipfman <- zm.ll(x, dist="Zipf-Man", s=2, b=3)

    ?zm.ll
    dist <- c(
        "Zeta Gamma 2",
        "Zeta",
        "Poisson",
        "Geometric",
        # "Zipf",
        # "Zipf-Mandelbrot",
        "Zeta Truncated"
    )

    logL <- c(
        2*minus_LL_zeta_gamma2(),
        attributes(summary(mle_zeta))$m2logL,
        attributes(summary(mle_poisson))$m2logL,
        attributes(summary(mle_geometric))$m2logL,
        # attributes(summary(mle_zipf))$m2logL,
        # attributes(summary(mle_zipfman))$m2logL,
        attributes(summary(mle_zeta_trunc))$m2logL
    )


    K <- c(
        0,
        1,
        1,
        1,
        2
    )

    coef <- c(
        list(""),
        list(attributes(summary(mle_zeta))$coef),
        list(attributes(summary(mle_poisson))$coef),
        list(attributes(summary(mle_geometric))$coef),
        # list(attributes(summary(mle_zipf))$coef),
        # list(attributes(summary(mle_zipfman))$coef),
        list(attributes(summary(mle_zeta_trunc))$coef)
    )

    #return(list(dist, logL, coef))

    #Calculate AIC
    aicV <- vector()
    for (i in 1:length(dist)) {
        aicV[i] <- get_AIC(logL[i],K[i],N)
    }

    #AICbest <- the one that minimizes AIC
    AICbest <- which.min(aicV)

    return(list(dist=dist[AICbest], params=coef[AICbest], candDist=dist, aic=aicV,deltaAic=aicV-aicV[AICbest]))
}



