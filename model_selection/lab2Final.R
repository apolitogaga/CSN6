# Cmplex Social Networks
# lab 2
# out degree Apolo Rosales Alberto Gutierrez
#install.packages("igraph")
library(igraph)
library(tolerance)
library(xtable)

write_summary <- function(language,file) {
  degree_sequence = read.table(file, header = FALSE)
  return (c(length(degree_sequence$V1),max(degree_sequence$V1),sum(degree_sequence$V1)/length(degree_sequence$V1),length(degree_sequence$V1)/sum(degree_sequence$V1) ))
}

source = read.table("list.txt",
                    header = TRUE,               # this is to indicate the first line of the file contains the names of the columns instead of the real data
                    as.is = c("language","file") # this is need to have the cells treated as real strings and not as categorial data.
)

names  = c("Arabic","Basque","Catalan","Chinese","Czech","English","Greek","Hungarian","Italian","Turkish")
shortnames  = c("ar","ba","ca","ch","cz","en","gr","hu","it","tu")
length(names)
languages = data.frame(10,col.names=names) #????
var=0
test= matrix(nrow=10,ncol=4)

#Read the data and do summary`
for (x in 1:nrow(source)) {
  temp <- write_summary(source$language[x], source$file[x])
  var[x] = read.table(source$file[x],header=FALSE)
  class(temp)
  for(y in 1:length(temp))
  {
    test[x,y] = temp[y]
  }
}

names(var) = shortnames


#Plot data
plotMe <-  function (list, name) {
  slist = sort(list, decreasing = TRUE)
  plot(slist, main = name, xlab = "log degree", ylab = "log Vertice #", log = "xy", type="l")
}
par(mfrow=(c(2,5)))
for(v in 1:length(var))
{
  plotMe(var[[v]], names[v])
}



### estimating parameters in R
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
  2 * sum(log(x)) + length(x) * log((pi^2)/6)
}


minus_LL_zeta_trunc <- function(gamma, Kmax) {
  H <- 0
  for(i in 1:Kmax) {
    H <- H + i^(-gamma)
  }
  gamma * sum(log(x)) + length(x) * log(H)
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
  mle_zeta_trunc <- mle(minus_LL_zeta_trunc,
                  start = list(gamma=gamma_t0,Kmax=length(x)),
                  method = "L-BFGS-B",
                  lower = c(1.000001))


  mle_zeta <- mle(minus_LL_zeta,
                  start = list(gamma=2),
                  method = "L-BFGS-B",
                  lower = c(1.000001))

  mle_poisson <-  mle(minus_LL_poisson,
                      start = list(lambda=mu),
                      method = "L-BFGS-B",
                      lower = c(0.000001))

  mle_geometric <- mle(minus_LL_geometric,
                       start = list(q=(1/mu)),
                       method = "L-BFGS-B",
                       lower = c(0),
                       upper = 0.999999)

  #mle_zipf <- mle(zipf_LL,
  #           start=list(s=1),
  #           method="L-BFGS-B",
  #           lower=c(1.0000001))
  mle_zipf <- zm.ll(x, dist="Zipf", s=1)

  #mle_zipfman <- mle(zipfMan_LL,
  #                start=list(s=1,q=1),
  #                method="L-BFGS-B",
  #                lower=c(1,0.0001))

  mle_zipfman <- zm.ll(x, dist="Zipf-Man", s=1, b=1)


  dist <- c(
    "Zeta Gamma 2",
    "Zeta",
    "Zeta Truncated",
    "Poisson",
    "Geometric",
    "Zipf",
    "Zipf-Mandelbrot"
  )

  logL <- c(
    2*minus_LL_zeta_gamma2(),
    attributes(summary(mle_zeta))$m2logL,
    attributes(summary(mle_zeta_trunc))$m2logL,
    attributes(summary(mle_poisson))$m2logL,
    attributes(summary(mle_geometric))$m2logL,
    attributes(summary(mle_zipf))$m2logL,
    attributes(summary(mle_zipfman))$m2logL
  )


  K <- c(
    0,
    1,
    2,
    1,
    1,
    1,
    2
  )

  coef <- c(
    list(""),
    list(attributes(summary(mle_zeta))$coef),
    list(attributes(summary(mle_zeta_trunc))$coef),
    list(attributes(summary(mle_poisson))$coef),
    list(attributes(summary(mle_geometric))$coef),
    list(attributes(summary(mle_zipf))$coef),
    list(attributes(summary(mle_zipfman))$coef)
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

#Testbench
x <- unlist(read.table("samples/data/sample_of_geometric_with_parameter_0.1.txt",header = FALSE))
findBestModel(3)

x <- unlist(read.table("./data/sample_of_geometric_with_parameter_0.2.txt",header = FALSE))
findBestModel(3)

x <- unlist(read.table("./data/sample_of_geometric_with_parameter_0.4.txt",header = FALSE))
findBestModel(3)

x <- unlist(read.table("./data/sample_of_geometric_with_parameter_0.8.txt",header = FALSE))
findBestModel(3)

#Error land.
x <- unlist(read.table("./data/sample_of_zeta_with_parameter_1.5.txt",header = FALSE))
findBestModel(3)
#f Error land.

x <- unlist(read.table("./data/sample_of_zeta_with_parameter_2.txt",header = FALSE))
findBestModel(3)

x <- unlist(read.table("./data/sample_of_zeta_with_parameter_2.5.txt",header = FALSE))
findBestModel(3)


x <- unlist(read.table("./data/sample_of_zeta_with_parameter_3.txt",header = FALSE))
findBestModel(3)



#Fit languages
gamma_t0 <- 2.6

mods <- sapply(var, function(x) {
  x <<- unlist(x)
  list(findBestModel(gamma_t0))
})

#View Distributions
(langDists <- sapply(mods, function(x) x$dist))

xtable(data.frame(language=names, distribution=langDists))
#View Parameters
langParams <- sapply(mods, function(x) x$params)


#Delta AIC table
delta <- t(sapply(mods, function(x) x$deltaAic))
rownames(delta) <- languages$col.names
colnames(delta) <- c(
                      "Zeta Gamma 2",
                      "Zeta",
                      "Zeta Truncated",
                      "Poisson",
                      "Geometric",
                      "Zipf",
                      "Zipf-Mandelbrot"
                    )

xtable(delta)


##PLOTS
## obtained from the parameter results table
param=matrix(nrow=10,ncol=2)
param = as.data.frame(param)
names(param)  = c("s","b")
row.names(param) = names
#ar
param[1,1]= 3.2707470#:06146179
param[1,2]= 3.2855440#:14222510
#basque
param[2,1]= 2.32434660#:05209573
param[2,2]= 0.84548560#:08970688
#CatalanZipf-Mandelbrot
param[3,1]= 2.5318790#:02817975
param[3,2]= 3.6176470#11182739
#ChineseZipf-Mandelbrot
param[4,1]= 2.3512140#:02351944
param[4,2]= 2.1030180#:07104917
#CzechZipf-Mandelbrot
param[5,1]= 2.7740590#:02543587
param[5,2]= 3.0865880#:07547503
#EnglishZipf-Mandelbrot
param[6,1]= 2.4147290#:02869153
param[6,2]= 3.8841360#:13133521
#GreekZipf-Mandelbrot
param[7,1]= 4.0661300#:1253984
param[7,2]= 6.6598760#:3741114
#HungarianZipf-Mandelbrot
param[8,1]= 2.9847160#:04050445
param[8,2]= 2.8957360#:09987030
#ItalianZipf-Mandelbrot
param[9,1]= 3.9887730#:1044263
param[9,2]= 6.3660660#:3091853
#TurkishZipf-Mandelbrot
param[10,1]= 4.1327480#:1056394
param[10,2]= 3.0691510#:1597280

####################################################################################
####################################################################################
####################################################################################
plotMe <-  function (list, name,s,b) {
  slist = sort(list, decreasing = TRUE)

  jpeg(paste("plots/",name,"fitted.jpg"))
  print(class(s))
  print(s)
  fitzim<- rzipfman(n = length(slist), s = s,b=b,N=length(slist))
  fitzim <- sort(fitzim,decreasing = TRUE)
  plot(slist, main = name, xlab = "log degree", ylab = "log Vertice #", log = "xy", type="l")
  lines(fitzim, col = 2, lwd = 2)
  dev.off()
}
