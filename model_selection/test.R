#
# funs = list(
#     list(f = nll_disp_poisson,  val = list(lambda=M/N),      low = c(0.0001),   up = c(Inf)),
#     list(f = nll_disp_geometric,val = list(q=N/M),           low = c(0.000001), up = c(0.999999)),
#     list(f = nll_zeta,          val = list(gamma=3),         low = c(1.0000001),up = c(Inf)),
#     list(f = nll_rt_zeta,       val = list(kmax=KM, gamma=3), low = c(N, 1.0), up = c(Inf, Inf)),
#     list(f = nll_altmann,       val = list(gamma=3, delta=0.001), low = c(0, 0), up = c(99, 99))
# )
# ?mle

language <-x1

lm(log(mean_length)~log(vertices), language)
with(language, leveneTest(mean_length,vertices,language))


###########################################################################
###############################   MODEL0   ################################
###########################################################################
nlMod0 <- function(language) {
    m0 <- as.formula("mean_length~a*(vertices)")

    #figure this out
    linear_model = lm(mean_length~vertices, language)
    ai = as.double(coef(linear_model)[2])

    #homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log((vertices/2)^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm = nls(m0,data=language, start = list(a = ai), trace = TRUE)
    param = c(coef(nonlm)["a"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
##############################   MODEL0+   ################################
###########################################################################
nlMod0p <- function(language) {
    m0 <- as.formula("mean_length~a*vertices+d")

    #figure this out
    linear_model = lm(mean_length~vertices, language)
    ai = as.double(coef(linear_model)[2])# vertices
    di = as.double(coef(linear_model)[1])# intercept
    #homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log((vertices/2)^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }
    nonlm = nls(m0,data=language, start = list(a = ai,d=di), trace = TRUE)
    param = c(coef(nonlm)["a"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

nlMod0p(language)

table(language$vertices)
###########################################################################
###############################   MODEL1   ################################
###########################################################################
nlMod1 <- function(language) {
    m1 <- as.formula("mean_length~a*(vertices)^(1/2)")

    #figure this out
    temp =log(language$vertices)/2.0
    linear_model = lm(log(mean_length)~temp, language)
    ai = exp(coef(linear_model))[2]

    ## homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log((vertices/2)^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm = nls(m1,data=language, start = list(a = ai), trace = TRUE)
    param = c(coef(nonlm)["a"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
###############################  MODEL1+   ################################
###########################################################################
nlMod1 <- function(language) {
    m1 <- as.formula("mean_length~a*(vertices)^(1/2)+d")

    #figure this out
    temp =log(language$vertices)/2.0
    linear_model = lm(log(mean_length)~temp, language)
    ai = exp(coef(linear_model))[2]
    di = 0

    ## homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log((vertices/2)^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm = nls(m1,data=language, start = list(a = ai,d=di), trace = TRUE)
    param = c(coef(nonlm)["a"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}
nlMod1(language)


###########################################################################
###############################   MODEL2   ################################
###########################################################################
nlMod2 <- function(language) {
    m2 <- as.formula("mean_length~a*vertices^b")

    linear_model = lm(log(mean_length)~log(vertices), language)
    ai = exp(coef(linear_model)[1])
    bi = coef(linear_model)[2]

    #homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*vertices^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm=nls(m2,data=language, start = list(a = ai, b = bi), trace = TRUE)
    param = c(coef(nonlm)["a"],coef(nonlm)["b"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
##############################   MODEL2+   ################################
###########################################################################
nlMod2p <- function(language) {
    m2 <- as.formula("mean_length~a*vertices^b+d")

    linear_model = lm(log(mean_length)~log(vertices), language)
    ai = exp(coef(linear_model)[1])
    bi = coef(linear_model)[2]
    d=0
    #homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*vertices^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm2=nls(m2,data=language, start = list(a = ai, b = bi,d=di), trace = TRUE,control=list(maxiter = 1000))
    param = c(coef(nonlm2)["a"],coef(nonlm2)["b"],coef(nonlm2)["d"])
    res = list("param" = param, "nonlm" = nonlm2)
    return(res)
}#3436.6903591    -0.0276991 -3162.1264523

###########################################################################
###############################   MODEL3   ################################
###########################################################################
nlMod3p <- function(language) {
    m3 <- as.formula("mean_length~a*exp(c*vertices)")

    linear_model = lm(log(mean_length)~vertices, language)
    ai = exp(coef(linear_model)[1])
    ci = coef(linear_model)[2]

    #homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*vertices^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm=nls(m3,data=language, start = list(a = ai, c = ci), trace = TRUE)
    param = c(coef(nonlm)["a"],coef(nonlm)["c"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
##############################   MODEL3+   ################################
###########################################################################
nlMod3 <- function(language) {
    m3 <- as.formula("mean_length~a*exp(c*vertices)+d")

    linear_model = lm(log(mean_length)~vertices, language)
    ai = exp(coef(linear_model)[1])
    ci = coef(linear_model)[2]
    di=2

    #homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*vertices^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm=nls(m3,data=language, start = list(a = ai, c = ci,d=di), trace = TRUE)
    param = c(coef(nonlm)["a"],coef(nonlm)["c"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

###########################################################################
###############################   MODEL4   ################################
###########################################################################
nlMod4 <- function(language) {
    m4 <- as.formula("mean_length~a*log(d+vertices)")

    linear_model = lm((exp(mean_length))~vertices, language)
    ai = 0
    di = coef(linear_model)[1]

    #homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*vertices^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm=nls(m4,data=language, start = list(a = ai, d = di), trace = TRUE)
    param = c(coef(nonlm)["a"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}
nlMod4(language)

###########################################################################
##############################   MODEL4+   ################################
###########################################################################
nlMod4p <- function(language) {
    m4 <- as.formula("mean_length~a*log(d1+vertices)+d2")

    linear_model = lm(exp(mean_length)~vertices, language)
    ai = 0
    di1 = coef(linear_model)[1]
    di2 = di1-1

    #homoscedasticity test if less than .05 we reject the null hypothesis
    # if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*vertices^b_initial))))$`Pr(>F)`[1] < 0.05){
    #     # language = aggLang(language)
    # }

    nonlm=nls(m4,data=language, start = list(a = ai, d1 = di1,d2=di2), trace = TRUE)
    param = c(coef(nonlm)["a"],coef(nonlm)["d"])
    res = list("param" = param, "nonlm" = nonlm)
    return(res)
}

nlMod4(language)



