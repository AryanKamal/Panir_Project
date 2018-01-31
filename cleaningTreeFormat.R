options(stringsAsFactors = F)

infoTree <- read.delim("2018 jan/Infomap Community Detection/new/germNet-FC1-subNetsNames.csv",sep = ",",header = T)



#for(i in 1:max(infoTree$V1)) {
#  n <- max(infoTree$V1[which(as.numeric(infoTree$V1) == i)])
#  for(j in 1:n) {
#    p <- max(infoTree$V2[which(as.numeric(infoTree$V2) == j)])
#    for(k in 1:p) {
#      v[length(v)+1] <- list(infoTree$names[which(as.numeric(infoTree$V3) == k)])
#    }
#  }
#}
v <- c()
name = c()
for(i in 1:max(infoTree $V1)) {
  n <- as.numeric(max(infoTree$V2[which(as.numeric(infoTree$V1) == i)]))
  for(j in 1:n) {
    new = list(infoTree$names[which(as.numeric(infoTree$V2) == j & as.numeric(infoTree$V1) == i)])
    name[length(name)+1] <- paste(as.character(i) , as.character(j) ,sep=":")
    v[length(v)+1] <- new
  }
}
names(v) <- name
max.length <- max(sapply(v, length))
## Add NA values to list elements
v <- lapply(v, function(x) { c(x, rep(NA, max.length-length(x)))})
## Rbind
do.call(cbind, v) -> t
colnames(t) <- names(v)
t= as.data.frame(t)


salam = "salam"