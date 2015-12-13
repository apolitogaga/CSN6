dba1 = as.numeric(unlist(read.table("data/BA_degree_1.txt",header = FALSE)))
dba10 = as.numeric(unlist(read.table("data/BA_degree_10.txt",header = FALSE)))
dba100 = as.numeric(unlist(read.table("data/BA_degree_100.txt",header = FALSE)))
dba1000 = as.numeric(unlist(read.table("data/BA_degree_1000.txt",header = FALSE)))
dbaFinal = as.numeric(unlist(read.table("data/BA_degree_sequence.txt",header = FALSE)))

dng1 = as.numeric(unlist(read.table("data/NG_degree_1.txt",header = FALSE)))
dng10 = as.numeric(unlist(read.table("data/NG_degree_10.txt",header = FALSE)))
dng100 = as.numeric(unlist(read.table("data/NG_degree_100.txt",header = FALSE)))
dng1000 = as.numeric(unlist(read.table("data/NG_degree_1000.txt",header = FALSE)))
dngFinal = as.numeric(unlist(read.table("data/NG_degree_sequence.txt",header = FALSE)))

x1=data.frame(mean_length= dngFinal,vertices=dng1)
x10= data.frame(mean_length= dngFinal,vertices=dng10)
x100= data.frame(mean_length= dngFinal,vertices=dng100)
x1000= data.frame(mean_length= dngFinal,vertices=dng1000)

languages=list(x1,x10,x100,x1000)

drg1 = as.numeric(unlist(read.table("data/RG_degree_1.txt",header = FALSE)))
drg10 = as.numeric(unlist(read.table("data/RG_degree_10.txt",header = FALSE)))
drg100 = as.numeric(unlist(read.table("data/RG_degree_100.txt",header = FALSE)))
drg1000 = as.numeric(unlist(read.table("data/RG_degree_1000.txt",header = FALSE)))
drgFinal = as.numeric(unlist(read.table("data/RG_degree_sequence.txt",header = FALSE)))

t2=read.table("data/test2.txt",header=F)
#we sort the data to be used easily
# d1 =  sort(d1,decreasing=T)
# d10 =  sort(d10,decreasing=T)
# d100 =  sort(d100,decreasing=T)
# d1000 =  sort(d1000,decreasing=T)

# we add them to a list to use the apply functions easily
dbaList  = list(dba1,dba10,dba100,dba1000,dbaFinal)
dngList  = list(dng1,dng10,dng100,dng1000,dngFinal)
drgList  = list(drg1,drg10,drg100,drg1000,drgFinal)



par(mfrow=c(2,2))
ly=c("d1","d10","d100","d1000","dF")
i=1
lapply(drgList,function(x){
    plot(x,ylab=ly[i])
    i<<-i+1
    }
)


lapply(dbaList,function(x){
    plot(x,ylab=ly[i])
    i<<-i+1
}
)

lapply(dngList,function(x){
    plot(x,ylab=ly[i])
    i<<-i+1
}
)




source("lab2LIB.r")
#Fit languages
gamma_t0 <- 1
mods <- sapply(dbaList, function(x) {
    x <<- unlist(x)
    list(findBestModel(gamma_t0))
})

#Delta AIC table
delta <- t(sapply(mods, function(x) x$deltaAic))
row.names(delta) <- names
colnames(delta) <- c(
    "Zeta Gamma 2",
    "Zeta",
    "Poisson",
    "Geometric",
    "Zeta Truncated"
)

print(xtable(delta, caption = "NG results"))

source("lab3LIB.r")

makeEverything(languages, names)

x<-dba1000
mle(minus_LL_zeta_trunc,
    start = list(Kmax=max(x),gamma=3),
    method = "L-BFGS-B",
    lower = c(N, 1.0000001),
    upper = c(Inf,Inf))

mle_zeta_trunc <- tryCatch(
    mle(minus_LL_zeta_trunc,
        start = list(Kmax=max(x),gamma=3),
        method = "L-BFGS-B",
        lower = c(N, 1.0000001),
        upper = c(Inf,Inf)),
    error=function(e) NULL)

#
# x<-dngFinal
# N=length(x)
# mu <- sum(x)/N
# mle_zeta_trunc <- tryCatch(
#     mle(minus_LL_zeta_trunc,
#         start = list(gamma=3,Kmax=N),
#         method = "L-BFGS-B",
#         lower = c(N, 1.0000001),
#         upper = c(Inf,Inf)),
#     error=function(e) NULL)
# mle_zeta_trunc

x<-dba1
N=length(x)
M_dash=sum(log(x))
mle(minus_log_like_zetaright,
    start = list(gamma=3,Kmax=max(x)),
    method = "L-BFGS-B",
    lower = c(1.0000001))

