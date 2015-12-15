library(car) # used for the homoscedasticity test
library(xtable)
library(lmtest)


aggLang <-function(lang){
    return(aggregate(lang, list(lang$time), mean))
}
# lm(log(degree)~log(time), network)
# with(network, leveneTest(degree,time,network))


###########################################################################
###############################   MODEL0   ################################
###########################################################################
nlMod0 <- function(network) {
    m0 <- as.formula("degree~a*(time)")
    linear_model = lm(degree~time, network)
    ai = as.double(coef(linear_model)[2])

    nonlm = nls(m0,data=network, start = list(a = ai))
    param = c(coef(nonlm)["a"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
##############################   MODEL0+   ################################
###########################################################################
nlMod0p <- function(network) {
    m0 <- as.formula("degree~a*time+d")

    #figure this out
    linear_model = lm(degree~time, network)
    ai = as.double(coef(linear_model)[2])# time
    di = as.double(coef(linear_model)[1])# intercept

    nonlm = nls(m0,data=network, start = list(a = ai,d=di))
    param = c(coef(nonlm)["a"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
###############################   MODEL1   ################################
###########################################################################
nlMod1 <- function(network) {
    m1 <- as.formula("degree~a*(time)^(1/2)")

    #figure this out
    temp =log(network$time)/2.0
    linear_model = lm(log(degree)~temp, network)
    ai = exp(coef(linear_model))[2]

    nonlm = nls(m1,data=network, start = list(a = ai))
    param = c(coef(nonlm)["a"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
###############################  MODEL1+   ################################
###########################################################################
nlMod1p <- function(network) {
    m1 <- as.formula("degree~a*(time)^(1/2)+d")

    #figure this out
    temp =log(network$time)/2.0
    linear_model = lm(log(degree)~temp, network)
    ai = exp(coef(linear_model))[2]
    di = 0

    nonlm = nls(m1,data=network, start = list(a = ai,d=di))
    param = c(coef(nonlm)["a"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}



###########################################################################
###############################   MODEL2   ################################
###########################################################################
nlMod2 <- function(network) {
    m2 <- as.formula("degree~a*time^b")

    linear_model = lm(log(degree)~log(time), network)
    ai = exp(coef(linear_model)[1])
    bi = coef(linear_model)[2]

    nonlm=nls(m2,data=network, start = list(a = ai, b = bi))
    param = c(coef(nonlm)["a"],coef(nonlm)["b"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
##############################   MODEL2+   ################################
###########################################################################
nlMod2p <- function(network){
    m2 <- as.formula("degree~a*time^b+d")

    linear_model = lm(log(degree)~log(time), network)
    ai = exp(coef(linear_model)[1])
    bi = coef(linear_model)[2]
    di=0

    nonlm2=nls(m2,data=network, start = list(a = ai, b = bi,d=di),control=list(maxiter = 2000))
    param = c(coef(nonlm2)["a"],coef(nonlm2)["b"],coef(nonlm2)["d"])
    res = list("param" = param, "nonlm" = nonlm2)
    return(res)
}

###########################################################################
###############################   MODEL3   ################################
###########################################################################
nlMod3 <- function(network) {
    m3 <- as.formula("degree~a*exp(c*time)")

    linear_model = lm(log(degree)~time, network)
    ai = exp(coef(linear_model)[1])
    ci = coef(linear_model)[2]

    nonlm=nls(m3,data=network, start = list(a = ai, c = ci))
    param = c(coef(nonlm)["a"],coef(nonlm)["c"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
##############################   MODEL3+   ################################
###########################################################################
nlMod3p <- function(network) {
    m3 <- as.formula("degree~a*exp(c*time)+d")

    linear_model = lm(log(degree)~time, network)
    ai = exp(coef(linear_model)[1])
    ci = coef(linear_model)[2]
    di=2

    nonlm=nls(m3,data=network, start = list(a = ai, c = ci,d=di))
    param = c(coef(nonlm)["a"],coef(nonlm)["c"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
###############################   MODEL4   ################################
###########################################################################
nlMod4 <- function(network) {
    m4 <- as.formula("degree~a*log(d+time)")

    linear_model = lm((exp(degree))~time, network)
    ai = 0
    di = coef(linear_model)[1]

    nonlm=nls(m4,data=network, start = list(a = ai, d = di))
    param = c(coef(nonlm)["a"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}


###########################################################################
##############################   MODEL4+   ################################
###########################################################################
nlMod4p <- function(network) {
    m4 <- as.formula("degree~a*log(d1+time)+d2")

    linear_model = lm(exp(degree)~time, network)
    ai = 0
    di1 = coef(linear_model)[1]
    di2 = di1-1

    nonlm=nls(m4,data=network, start = list(a = ai, d1 = di1,d2=di2))
    param = c(coef(nonlm)["a"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}


###########################################################################
makeEverything <- function(lang, langnames, capName=""){
    nls_lang_m0 = sapply(lang, nlMod0)
    param_m0 = as.data.frame(nls_lang_m0[1,])
    nls_m0 = nls_lang_m0[2,]
    nls_lang_m0p = sapply(lang, nlMod0p)
    param_m0p = as.data.frame(nls_lang_m0p [1,])
    nls_m0p = nls_lang_m0p[2,]
    nls_lang_m1 = sapply(lang, nlMod1)
    param_m1 = as.data.frame(nls_lang_m1 [1,])
    nls_m1 = nls_lang_m1[2,]
    nls_lang_m2 = sapply(lang, nlMod2)
    param_m2 = as.data.frame(nls_lang_m2 [1,])
    nls_m2 = nls_lang_m2[2,]
    nls_lang_m3 = sapply(lang, nlMod3)
    param_m3 = as.data.frame(nls_lang_m3 [1,])
    nls_m3 = nls_lang_m3[2,]
    nls_lang_m1p = sapply(lang, nlMod1p)
    param_m1p = as.data.frame(nls_lang_m1p [1,])
    nls_m1p = nls_lang_m1p[2,]

    # nls_lang_m2p = sapply(lang, nlMod2p)
    # param_m2p = as.data.frame(nls_lang_m2p [1,])
    # nls_m2p = nls_lang_m2p[2,]
    # nls_lang_m3p = sapply(lang[2], nlMod3p)

    listList = list(nls_m0,nls_m1,nls_m2,nls_m3,nls_m0p,nls_m1p)
    # listList = list(nls_m0,nls_m0p)
    table2 =  as.data.frame(matrix(ncol=10,nrow=5))
    # View(table2)
    row.names(table2) <- langnames
    colnames(table2) <-c("0a","0pa","0pb","1a", "1pb",  "1pd", "2a",  "2b", "3a",  "3c" )

    # add 0 and 0p
    table2[,1]=t(param_m0)
    table2[,2]=t(param_m0p[1,])
    table2[,3]=t(param_m0p[2,])

    # # add 1 and 1p
    table2[,4]=t(param_m1)
    table2[,5]=t(param_m1p[1,])
    table2[,6]=t(param_m1p[2,])
    #
    # # add 2 and 2p
    table2[,7]=t(param_m2[1,])
    table2[,8]=t(param_m2[2,])
    # table2[,6]=t(param_m2p[1,])
    # table2[,7]=t(param_m2p[2,])
    # table2[,8]=t(param_m2p[3,])
    # # add 3
    table2[,9]=t(param_m3[1,])
    table2[,10]=t(param_m3[2,])
    #
    print("################ Best parameters table ##################")
    print(xtable(table2,digits = 5, caption=paste("Best parameters",capName) ))
    #
    # Calculating model 0 values
    model0 <- function(language){
        RSS <- sum((language$mean_length-(language$vertices+1)/3)^2)
        n <- length(language$vertices)
        p <- 0
        s <- sqrt(RSS/(n - p))
        AIC <- n*log(2*pi) + n*log(RSS/n) + n + 2*(p + 1)
        return(list("s"=s, "AIC"=AIC))
    }
    model0values = sapply(lang, model0)
    #
    # S Table
    svalues_m0 <- sapply(nls_m0, function(nls){sqrt(deviance(nls)/df.residual(nls))})
    svalues_m1 <- sapply(nls_m1, function(nls){sqrt(deviance(nls)/df.residual(nls))})
    svalues_m2 <- sapply(nls_m2, function(nls){sqrt(deviance(nls)/df.residual(nls))})
    svalues_m3 <- sapply(nls_m3, function(nls){sqrt(deviance(nls)/df.residual(nls))})
    svalues_m0p <- sapply(nls_m0p, function(nls){sqrt(deviance(nls)/df.residual(nls))})
    svalues_m1p <- sapply(nls_m1p, function(nls){sqrt(deviance(nls)/df.residual(nls))})
    # svalues_m2p <- sapply(nls_m2p, function(nls){sqrt(deviance(nls)/df.residual(nls))})

    tableS =  as.data.frame(matrix(ncol=6,nrow=5))
    row.names(tableS) <- langnames
    colnames(tableS) <-c("0", "1", "2", "3","0p","1p")

    # tableS[,1]=unlist(model0values[1,])
    tableS[,1]=svalues_m0
    tableS[,2]=svalues_m1
    tableS[,3]=svalues_m2
    tableS[,4]=svalues_m3
    tableS[,5]=svalues_m0p
    tableS[,6]=svalues_m1p
    # tableS[,7]=svalues_m2p

    # tableS
    # print("################ S TABLE ##################")
    print(xtable(tableS,caption=paste("S table ",capName)))
    #
    # # AIC Table
    aic_m0 <- sapply(nls_m0, function(nls){AIC(nls)})
    aic_m1 <- sapply(nls_m1, function(nls){AIC(nls)})
    aic_m2 <- sapply(nls_m2, function(nls){AIC(nls)})
    aic_m3 <- sapply(nls_m3, function(nls){AIC(nls)})
    aic_m0p <- sapply(nls_m0p, function(nls){AIC(nls)})
    aic_m1p <- sapply(nls_m1p, function(nls){AIC(nls)})
    # aic_m2p <- sapply(nls_m2p, function(nls){AIC(nls)})
    #
    tableAIC =  as.data.frame(matrix(ncol=6,nrow=5))
    row.names(tableAIC) <- langnames
    colnames(tableAIC) <-c("0", "1", "2", "3","0p", "1p")
    #
    tableAIC[,1]=aic_m0
    tableAIC[,2]=aic_m1
    tableAIC[,3]=aic_m2
    tableAIC[,4]=aic_m3
    tableAIC[,5]=aic_m0p
    tableAIC[,6]=aic_m1p
    # tableAIC[,7]=aic_m2p
    #
    # print("################ AIC TABLE ##################")
    print(xtable(tableAIC, caption=paste("AIC",capName)))
    #
    # # Delta AIC Table
    tableDeltaAIC = tableAIC
    tableDeltaAIC[1,] = abs(tableDeltaAIC[1,] - min(tableDeltaAIC[1,]))
    tableDeltaAIC[2,] = abs(tableDeltaAIC[2,] - min(tableDeltaAIC[2,]))
    tableDeltaAIC[3,] = abs(tableDeltaAIC[3,] - min(tableDeltaAIC[3,]))
    tableDeltaAIC[4,] = abs(tableDeltaAIC[4,] - min(tableDeltaAIC[4,]))
    tableDeltaAIC[5,] = abs(tableDeltaAIC[5,] - min(tableDeltaAIC[5,]))
    # # tableDeltaAIC[6,] = abs(tableDeltaAIC[6,] - min(tableDeltaAIC[6,]))
    # # tableDeltaAIC[7,] = abs(tableDeltaAIC[7,] - min(tableDeltaAIC[7,]))
    # # tableDeltaAIC[8,] = abs(tableDeltaAIC[8,] - min(tableDeltaAIC[8,]))
    # # tableDeltaAIC[9,] = abs(tableDeltaAIC[9,] - min(tableDeltaAIC[9,]))
    # # tableDeltaAIC[10,] = abs(tableDeltaAIC[10,] - min(tableDeltaAIC[10,]))
    #
    # print("################ DELTA AIC TABLE ##################")
    print(xtable(tableDeltaAIC,caption=paste("Delta AIC ",capName)))
    #
    #Final Plotting with best fits

    for(x in 1:length(lang)){
        minColumn = which.min(tableDeltaAIC[x,])
        if(minColumn == 2) {
            nonlinear_model = nls_m1[[x]]
        }
        if(minColumn == 3) {
            nonlinear_model = nls_m2[[x]]
        }
        if(minColumn == 4) {
            nonlinear_model = nls_m3[[x]]
        }
        if(minColumn == 5) {
            nonlinear_model = nls_m0p[[x]]
        }
        if(minColumn == 6) {
            nonlinear_model = nls_m1p[[x]]
        }
        # if(minColumn == 6) {
        #     nonlinear_model = nls_m2p[[x]]
        # }

        if(length(log(fitted(nonlinear_model))) < length(log(lang[[x]]$time)))
        {
            lang[[x]] = aggLang(lang[[x]])
        }
        p4Name = paste("plots/", gsub(" ","",capName),"FinalFIT",langnames[x],".jpg",sep = "")
        jpeg(p4Name)
        plot(log(lang[[x]]$time), log(lang[[x]]$degree), main = paste(" Land for",capName ,langnames[x]),
             xlab = "log(time)", ylab = "log(degree)")
        lines(log(lang[[x]]$time), log(fitted(nonlinear_model)), col = "green")
        dev.off()
    }

}

