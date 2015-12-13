#libraries
# install.packages("Rcmdr")
library(car) # used for the homoscedasticity test
library(xtable)
library(lmtest)
### Read the files from the data folders
#setwd("~/Dropbox/Homework 3")
ar = read.table("./data/Arabic_dependency_tree_metrics.txt", header = FALSE)
ba = read.table("./data/Basque_dependency_tree_metrics.txt", header = FALSE)
cat = read.table("./data/Catalan_dependency_tree_metrics.txt", header = FALSE)
ch = read.table("./data/Chinese_dependency_tree_metrics.txt", header = FALSE)
cz = read.table("./data/Czech_dependency_tree_metrics.txt", header = FALSE)
en = read.table("./data/English_dependency_tree_metrics.txt", header = FALSE)
gr = read.table("./data/Greek_dependency_tree_metrics.txt", header = FALSE)
hu = read.table("./data/Hungarian_dependency_tree_metrics.txt", header = FALSE)
it = read.table("./data/Italian_dependency_tree_metrics.txt", header = FALSE)
tu = read.table("./data/Turkish_dependency_tree_metrics.txt", header = FALSE)

#setting column names for easy access
colnames(ar) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(ba) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(cat) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(ch) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(cz) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(en) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(gr) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(hu) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(it) <- c("vertices","degree_2nd_moment", "mean_length")
colnames(tu) <- c("vertices","degree_2nd_moment", "mean_length")
# ar = ar[-2]
# ba = ba[-2]
# cat = cat[-2]
# ch = ch[-2]
# cz = cz[-2]
# en = en[-2]
# gr = gr[-2]
# hu = hu[-2]
# it = it[-2]
# tu = tu[-2]

#putting all languages in a list to access through

lang = list(ar,ba,cat,ch,cz,en,gr,hu,it,tu)

langnames  = c(
  "Arabic","Basque","Catalan","Chinese","Czech","English","Greek","Hungarian","Italian","Turkish"
)
langshortnames  = c("ar","ba","ca","ch","cz","en","gr","hu","it","tu")

names(lang) <- langshortnames

#Function for checking if the data conforms with the necessary conditions
validity <- function(language) {
  n = nrow(language[1])
  k = language[,2]
  d = language[,3]
  
  for (x in 1:n) {
    if (n / (8 * (n - 1)) * k[x] + (1 / 2) >= d[x]) {
      print(paste("Error1: ", abs (n / (8 * (
        n - 1
      )) * k[x] + (1 / 2) - d[x]), "on value ", x))
    }
    
    if (d[x] <= n - 1) {
      
    } else {
      print(paste("Error2: ", abs (n - 1 - d[x]), "on value ", x))
    }
  }
  
}

#plotting a language's mean dependency length in relation to a word's number of vertices in normal and logarithmic scale
plotFirst <- function(toPlot, name, extra=""){
  plotAgg
  p1Name =  paste("plots/",name,".jpg",sep="")
  p2Name = paste("plots/log",name,".jpg",sep="")
  jpeg(p1Name)
  plot(toPlot$vertices, toPlot$mean_length, main=name,
       xlab = "vertices", ylab = paste(extra,"mean dependency length",sep=""))
  dev.off()
  jpeg(p2Name)
  plot(log(toPlot$vertices), log(toPlot$mean_length), main=paste("log",name,sep=""),
       xlab = "log(vertices)", ylab = paste("log(",extra,"mean dependency length)",sep=""))
  dev.off()
}
aggLang <-function(lang){
  return(aggregate(lang, list(lang$vertices), mean))
}

plotAgg<-function(lang,name){
  p3Name = paste("plots/realScaling",name,".jpg",sep = "")
  mean_lang = aggLang(lang)
  plotFirst(mean_lang, paste("aggregate",name), "mean")
  jpeg(p3Name)
  plot(log(lang$vertices), log(lang$mean_length), main=paste(name,"scaling",sep=''),
       xlab = "vertices", ylab = "mean dependency length")
  lines(log(mean_lang$vertices),log(mean_lang$mean_length), col = "green")
  lines(log(mean_lang$vertices),log((mean_lang$vertices+1)/3), col = "red")
  dev.off()
}

# plot the real scaling of D and the linear arrangement.
# plotAll <- function(language, name){
#   plotFirst(language, name)
#   plotAgg(language, name)
# }

# #preliminary plots for all languages
# for(i in 1:length(lang)){
#   plotAll(lang[[i]],langnames[i])
# }

#function for plotting table 1
write_table1 <- function(languageName,language) {
  N <- length(language[[1]][[1]])
  un <- mean(language[[1]][[1]])
  on <- sd(language[[1]][[1]])
  ud <- mean(language[[1]][[3]])
  od <- sd(language[[1]][[3]])
  return (c( N, un, on, ud, od))
}

table1 = as.data.frame(matrix(nrow=10,ncol=5))
row.names(table1) <- langnames
colnames(table1)<-c("N","un","on","ud", "od")

for (i in 1:length(lang)) {
 table1[i,]<- write_table1(langnames[i], lang[i])
  
}
# print table 1
print("################ TABLE 1 ##################")
print(xtable(table1))

#model 0  -> f(n)=n/3+1/3       double log: b*log(n)-b*log(2)
#model 1  -> f(n)=(n/2)^b       double log: 
#model 2  -> f(n)=an^b          double log: 
#model 3  -> f(n)=ae^(cn)       double log: 
#model 1+ -> f(n)=(n/2)^b+d     double log: 
#model 2+ -> f(n)=an^b+d        double log: 

#EXTRA WORK
#model 3+ -> f(n)=ae^(cn)+d 
#model 4  -> f(n)=an^b*e^(cn)
#model 4+ -> f(n)=an^b*e^(cn)+d

###########################################################################
# MODEL 1

model1 <- function(language) {
  m1 <- as.formula("mean_length~(vertices/2)^b")
  
  #figure this out
  linear_model = lm(log(mean_length)~log(vertices), language)
  b_initial = coef(linear_model)[1]
  
  #homoscedasticity test if less than .05 we reject the null hypothesis
  if(with(language, leveneTest(log(mean_length), as.factor(log((vertices/2)^b_initial))))$`Pr(>F)`[1] < 0.05){
    # language = aggLang(language)
  }
  
  nonlm = nls(m1,data=language, start = list(b = b_initial), trace = TRUE)
  param = c(coef(nonlm)["b"])
  res = list("param" = param, "nonlm" = nonlm)
  return(res)
}
###########################################################################
# MODEL 2
model2 <- function(language) {
  m2 <- as.formula("mean_length~a*vertices^b")
  
  linear_model = lm(log(degree_2nd_moment)~log(vertices), language)
  a_initial = exp(coef(linear_model)[1])
  b_initial = coef(linear_model)[2]
  
  #homoscedasticity test if less than .05 we reject the null hypothesis
  if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*vertices^b_initial))))$`Pr(>F)`[1] < 0.05){
    # language = aggLang(language)
  }
  
  nonlm=nls(m2,data=language, start = list(a = a_initial, b = b_initial), trace = TRUE)
  param = c(coef(nonlm)["a"],coef(nonlm)["b"])
  res = list("param" = param, "nonlm" = nonlm)
  return(res)
}
###########################################################################
# MODEL 3
model3 <- function(language) {
  n <- nrow(language)
  m3 <- as.formula("mean_length~a*exp(c*vertices)")
  
  linear_model = lm(log(mean_length)~vertices, language)
  a_initial = exp(coef(linear_model)[1])
  c_initial = coef(linear_model)[2]
  
  #homoscedasticity test if less than .05 we reject the null hypothesis
  if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*exp(c_initial*vertices)))))$`Pr(>F)`[1] < 0.05){
    # language = aggLang(language)
  }
  
  nonlm=nls(m3,data=language, start = list(a = a_initial, c = c_initial), trace = TRUE)
  param = c(coef(nonlm)["a"],coef(nonlm)["c"])
  res = list("param" = param, "nonlm" = nonlm)
  return(res)
}
# 
###########################################################################
# # MODEL 4
# model4 <- function(language) {
#   m4 <- as.formula("mean_length~a*n^b*e^(c*n)")
#   
#   linear_model = lm(log(mean_length)~vertices, language)
#   a_initial = exp(coef(linear_model)[1])
#   c_initial = coef(linear_model)[2]
#   
#   nls(m3,data=language, start = list(a = a_initial, c = c_initial), trace = TRUE)
# }
log(ar$mean_length)
log(ar$vertices)
###########################################################################
# MODEL 1+
model1p <- function(language)
{
  n <- nrow(language)
  m1p <- as.formula("mean_length~(vertices/2)^b+d")
  
  linear_model = lm(log(mean_length)~log(vertices), data=language)
  b_initial = coef(linear_model)[2]
  d_initial = 1
  
  #homoscedasticity test if less than .05 we reject the null hypothesis
  if(with(language, leveneTest(log(mean_length), as.factor(log((vertices/2)^b_initial+d_initial))))$`Pr(>F)`[1] < 0.05){
    # language = aggLang(language)
  }
  
  nonlm = nls(m1p,data=language, start = list(b = b_initial, d = d_initial), trace = TRUE)
  param = c(coef(nonlm)["b"],coef(nonlm)["d"])
  res = list("param" = param, "nonlm" = nonlm)
  return(res)
}
###########################################################################
# MODEL 2+

model2p <- function(language) {
  n <- nrow(language)
  m2p <- as.formula("mean_length~a*vertices^b+d")
  
  linear_model = lm(log(mean_length)~log(vertices), language)
  a_initial = exp(coef(linear_model)[1])
  b_initial = coef(linear_model)[2]
  d_initial = 0
  
  #homoscedasticity test if less than .05 we reject the null hypothesis
  if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*vertices^b_initial+d_initial))))$`Pr(>F)`[1] < 0.05){
    # language = aggLang(language)
  }
  
  nonlm = nls(m2p, data=language, start = list(a = a_initial, b = b_initial, d = d_initial), trace = TRUE)
  param =c(coef(nonlm)["a"],coef(nonlm)["b"],coef(nonlm)["d"])
  res = list("param" = param, "nonlm" = nonlm)
  return(res)
}
###########################################################################
# MODEL 3+
model3p <- function(language) {
  n <- nrow(language)
  m3p <- as.formula("mean_length~a*exp(c*vertices) + d ")
  
  linear_model = lm(log(mean_length)~log(vertices), language)
  a_initial = exp(coef(linear_model)[1])
  c_initial = coef(linear_model)[2]
  d_initial = 2
  
  #homoscedasticity test if less than .05 we reject the null hypothesis
  if(with(language, leveneTest(log(mean_length), as.factor(log(a_initial*exp(c_initial*vertices)))))$`Pr(>F)`[1] < 0.05){
    # language = aggLang(language)
  }
  
  nonlm = nls(m3p, data=language, start = list(a = a_initial, c = c_initial, d = d_initial), trace = TRUE)
  param = c(coef(nonlm)["a"],coef(nonlm)["c"],coef(nonlm)["d"])
  res = list("param" = param, "nonlm" = nonlm)
  return(res)
}

###########################################################################
makeEverything <- function(lang, langnames){
  nls_lang_m1 = sapply(lang, model1)
  param_m1 = as.data.frame(nls_lang_m1 [1,])
  nls_m1 = nls_lang_m1[2,]
  nls_lang_m2 = sapply(lang, model2)
  param_m2 = as.data.frame(nls_lang_m2 [1,])
  nls_m2 = nls_lang_m2[2,]
  nls_lang_m3 = sapply(lang, model3)
  param_m3 = as.data.frame(nls_lang_m3 [1,])
  nls_m3 = nls_lang_m3[2,]
  nls_lang_m1p = sapply(lang, model1p)
  param_m1p = as.data.frame(nls_lang_m1p [1,])
  nls_m1p = nls_lang_m1p[2,]
  nls_lang_m2p = sapply(lang, model2p)
  param_m2p = as.data.frame(nls_lang_m2p [1,])
  nls_m2p = nls_lang_m2p[2,]
  #nls_lang_m3p = sapply(lang[2], model3p)
  
  listList = list(nls_m1,nls_m2,nls_m3,nls_m1p,nls_m2p)
  table2 =  as.data.frame(matrix(ncol=10,nrow=10))
  row.names(table2) <- langnames
  colnames(table2) <-c("1b",  "1pb",  "1pd", "2a",  "2b",  "2pa",  "2pb",  "2pd",  "3a",  "3c" )
  # add 1 and 1p
  table2[,1]=t(param_m1)
  table2[,2]=t(param_m1p[1,])
  table2[,3]=t(param_m1p[2,])
  
  # add 2 and 2p
  table2[,4]=t(param_m2[1,])
  table2[,5]=t(param_m2[2,])
  table2[,6]=t(param_m2p[1,])
  table2[,7]=t(param_m2p[2,])
  table2[,8]=t(param_m2p[3,])
  # add 3
  table2[,9]=t(param_m3[1,])
  table2[,10]=t(param_m3[2,])
  
  print("################ Best parameters table ##################")
  print(xtable(table2))
  
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
  
  # S Table
  svalues_m1 <- sapply(nls_m1, function(nls){sqrt(deviance(nls)/df.residual(nls))})
  svalues_m2 <- sapply(nls_m2, function(nls){sqrt(deviance(nls)/df.residual(nls))})
  svalues_m3 <- sapply(nls_m3, function(nls){sqrt(deviance(nls)/df.residual(nls))})
  svalues_m1p <- sapply(nls_m1p, function(nls){sqrt(deviance(nls)/df.residual(nls))})
  svalues_m2p <- sapply(nls_m2p, function(nls){sqrt(deviance(nls)/df.residual(nls))})
  
  tableS =  as.data.frame(matrix(ncol=6,nrow=10))
  row.names(tableS) <- langnames
  colnames(tableS) <-c("0", "1", "2", "3", "1p", "2p")
  
  tableS[,1]=unlist(model0values[1,])
  tableS[,2]=svalues_m1
  tableS[,3]=svalues_m2
  tableS[,4]=svalues_m3
  tableS[,5]=svalues_m1p
  tableS[,6]=svalues_m2p
  
  tableS
  print("################ S TABLE ##################")
  print(xtable(tableS))
  
  # AIC Table
  aic_m1 <- sapply(nls_m1, function(nls){AIC(nls)})
  aic_m2 <- sapply(nls_m2, function(nls){AIC(nls)})
  aic_m3 <- sapply(nls_m3, function(nls){AIC(nls)})
  aic_m1p <- sapply(nls_m1p, function(nls){AIC(nls)})
  aic_m2p <- sapply(nls_m2p, function(nls){AIC(nls)})
  
  tableAIC =  as.data.frame(matrix(ncol=6,nrow=10))
  row.names(tableAIC) <- langnames
  colnames(tableAIC) <-c("0", "1", "2", "3", "1p", "2p")
  
  tableAIC[,1]=unlist(model0values[2,])
  tableAIC[,2]=aic_m1
  tableAIC[,3]=aic_m2
  tableAIC[,4]=aic_m3
  tableAIC[,5]=aic_m1p
  tableAIC[,6]=aic_m2p
  
  print("################ AIC TABLE ##################")
  print(xtable(tableAIC))
  
  # Delta AIC Table
  tableDeltaAIC = tableAIC
  tableDeltaAIC[1,] = abs(tableDeltaAIC[1,] - min(tableDeltaAIC[1,]))
  tableDeltaAIC[2,] = abs(tableDeltaAIC[2,] - min(tableDeltaAIC[2,]))
  tableDeltaAIC[3,] = abs(tableDeltaAIC[3,] - min(tableDeltaAIC[3,]))
  tableDeltaAIC[4,] = abs(tableDeltaAIC[4,] - min(tableDeltaAIC[4,]))
  tableDeltaAIC[5,] = abs(tableDeltaAIC[5,] - min(tableDeltaAIC[5,]))
  tableDeltaAIC[6,] = abs(tableDeltaAIC[6,] - min(tableDeltaAIC[6,]))
  tableDeltaAIC[7,] = abs(tableDeltaAIC[7,] - min(tableDeltaAIC[7,]))
  tableDeltaAIC[8,] = abs(tableDeltaAIC[8,] - min(tableDeltaAIC[8,]))
  tableDeltaAIC[9,] = abs(tableDeltaAIC[9,] - min(tableDeltaAIC[9,]))
  tableDeltaAIC[10,] = abs(tableDeltaAIC[10,] - min(tableDeltaAIC[10,]))
  
  print("################ DELTA AIC TABLE ##################")
  print(xtable(tableDeltaAIC))
   
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
      nonlinear_model = nls_m1p[[x]]
    }
    if(minColumn == 6) {
      nonlinear_model = nls_m2p[[x]]
    }
    
    if(length(log(fitted(nonlinear_model))) < length(log(lang[[x]]$vertices)))
    {
      lang[[x]] = aggLang(lang[[x]])
    }
    p4Name = paste("plots/finalFIT",langnames[x],".jpg",sep = "")
    jpeg(p4Name)
    plot(log(lang[[x]]$vertices), log(lang[[x]]$mean_length), main = langnames[x],
         xlab = "log(vertices)", ylab = "log(mean dependency length)")
    lines(log(lang[[x]]$vertices), log(fitted(nonlinear_model)), col = "green")
    dev.off()
  }
  
}
makeEverything(lang, langnames)
