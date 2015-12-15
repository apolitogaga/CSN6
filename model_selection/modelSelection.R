dba1 = as.numeric(unlist(read.table("data/BA_degree_1.txt",header = FALSE)))
dba1 = cbind(dba1,seq(1,length(dba1)))
colnames(dba1) <-c("degree","time")
dba1=as.data.frame(dba1)
dba10 = as.numeric(unlist(read.table("data/BA_degree_10.txt",header = FALSE)))
dba10 = cbind(dba10,seq(1,length(dba10)))
colnames(dba10) <-c("degree","time")
dba10=as.data.frame(dba10)
dba100 = as.numeric(unlist(read.table("data/BA_degree_100.txt",header = FALSE)))
dba100 = cbind(dba100,seq(1,length(dba100)))
colnames(dba100) <-c("degree","time")
dba100=as.data.frame(dba100)
dba1000 = as.numeric(unlist(read.table("data/BA_degree_1000.txt",header = FALSE)))
dba1000 = cbind(dba1000,seq(1,length(dba1000)))
colnames(dba1000) <-c("degree","time")
dba1000=as.data.frame(dba1000)
dbaFinal = as.numeric(unlist(read.table("data/BA_degree_sequence.txt",header = FALSE)))
dbaFinal = cbind(dbaFinal,seq(1,length(dbaFinal)))
colnames(dbaFinal) <-c("degree","time")
dbaFinal = as.data.frame(dbaFinal)


dng1 = as.numeric(unlist(read.table("data/NG_degree_1.txt",header = FALSE)))
dng1 = cbind(dng1,seq(1,length(dng1)))
colnames(dng1) <-c("degree","time")
dng1=as.data.frame(dng1)
dng10 = as.numeric(unlist(read.table("data/NG_degree_10.txt",header = FALSE)))
dng10 = cbind(dng10,seq(1,length(dng10)))
colnames(dng10) <-c("degree","time")
dng10=as.data.frame(dng10)
dng100 = as.numeric(unlist(read.table("data/NG_degree_100.txt",header = FALSE)))
dng100 = cbind(dng100,seq(1,length(dng100)))
colnames(dng100) <-c("degree","time")
dng1000 = as.numeric(unlist(read.table("data/NG_degree_1000.txt",header = FALSE)))
dng1000 = cbind(dng1000,seq(1,length(dng1000)))
colnames(dng1000) <-c("degree","time")
dngFinal = as.numeric(unlist(read.table("data/NG_degree_sequence.txt",header = FALSE)))
dngFinal = cbind(dngFinal,seq(1,length(dngFinal)))
colnames(dngFinal) <-c("degree","time")
dng100=as.data.frame(dng100)
dng1000=as.data.frame(dng1000)
dngFinal=as.data.frame(dngFinal)

drg1 = as.numeric(unlist(read.table("data/RG_degree_1.txt",header = FALSE)))
drg1 = cbind(drg1,seq(1,length(drg1)))
drg10 = as.numeric(unlist(read.table("data/RG_degree_10.txt",header = FALSE)))
drg10 = cbind(drg10,seq(1,length(drg10)))
drg100 = as.numeric(unlist(read.table("data/RG_degree_100.txt",header = FALSE)))
drg100 = cbind(drg100,seq(1,length(drg100)))
drg1000 = as.numeric(unlist(read.table("data/RG_degree_1000.txt",header = FALSE)))
drg1000 = cbind(drg1000,seq(1,length(drg1000)))
drgFinal = as.numeric(unlist(read.table("data/RG_degree_sequence.txt",header = FALSE)))
drgFinal = cbind(drgFinal,seq(1,length(drgFinal)))

colnames(drg1) <-c("degree","time")
colnames(drg10) <-c("degree","time")
colnames(drg100) <-c("degree","time")
colnames(drg1000) <-c("degree","time")
colnames(drgFinal) <-c("degree","time")

drg1=as.data.frame(drg1)
drg10=as.data.frame(drg10)
drg100=as.data.frame(drg100)
drg1000=as.data.frame(drg1000)
drgFinal=as.data.frame(drgFinal)

# we add them to a list to use the apply functions easily
dbaList  = list(dba1,dba10,dba100,dba1000,dbaFinal)
dngList  = list(dng1,dng10,dng100,dng1000,dngFinal)
drgList  = list(drg1,drg10,drg100,drg1000,drgFinal)
names=c("d1","d10","d100","d1000","dF")



curve(dba1$degree[1] * sqrt(x / 1) * sqrt(1),from = 1, to = dba1$time[length(dba1$time)], col = "green")
curve(dba1$degree[1] * sqrt(x / 10) * sqrt(10),from = 10, to = dba1$time[length(dba1$time)], col = "red", add = TRUE)
curve(dba1$degree[1] * sqrt(x / 100) * sqrt(100),from = 100, to = dba1$time[length(dba1$time)], col = "blue", add = TRUE)
curve(dba1$degree[1] * sqrt(x / 1000) * sqrt(1000),from = 1000, to = dba1$time[length(dba1$time)], col = "yellow", add = TRUE)
curve(dba1$degree[1] * sqrt(x),from = 1, to = dba1$time[length(dba1$time)], col = "red", type = "l", add = TRUE)

par(mfrow=c(2,2))
ly=c("d1","d10","d100","d1000","dF")
i=1



source("lab2LIB.r")
test<-function(list,name=""){

#Fit languages
gamma_t0 <- 1
mods <- sapply(list, function(x) {
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
print(xtable(delta, caption = paste(name,"results")))
}

test(dngList,"No Growth")
test(drgList,"Random Growth")
test(dbaList,"No Growth")



source("lab3LIB.R")
makeEverything(dngList, names,"No Growth")
makeEverything(drgList, names,"Random Growth")
makeEverything(dbaList, names,"Barabasi-Albert")


