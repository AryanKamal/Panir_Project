options(stringsAsFactors = F)
library(igraph)

DElist1 <- read.delim("afterDefence/2data-cb-easyDEex.txt",sep = "\t",header = T)
colnames(DElist1)[1] <- "Name"
DElist1 <- DElist1[,1]

DElist2 <- read.delim("afterDefence/2data-cb-FC2DEex.txt",sep = "\t",header = T)
colnames(DElist2)[1] <- "Name"
DElist2 <- DElist2[,1]


cellnet <- read.delim("C:/Users/Asus/Desktop/Repeat-Jun2017/CellNet/Mouse_Big_GRN.txt",
                      sep = " ",header = T)
cellnet <- cellnet[,c(2,1,3:6)]

cellnetGraph<-graph_from_data_frame(cellnet,directed=T)


CommonList1<-intersect(unlist(DElist1),V(cellnetGraph)$name)
CommonList2<-intersect(unlist(DElist2),V(cellnetGraph)$name)


indGraph1<-induced_subgraph(cellnetGraph,unlist(CommonList1))
indGraph2<-induced_subgraph(cellnetGraph,unlist(CommonList2))

subNet1<-as_long_data_frame(indGraph1)  
subNet2<-as_long_data_frame(indGraph2)  


germ1 <- subNet1[which(subNet1$type == "germ_esc" | subNet1$type == "all"),]
germ2 <- subNet2[which(subNet2$type == "germ_esc" | subNet2$type == "all"),]

all1link <- subNet1[,c(1,2,4)]
all1link$corr <- abs(all1link$corr)

all2link <- subNet2[,c(1,2,4)]
all2link$corr <- abs(all2link$corr)

germ1link <- germ1[,c(1,2,4)]
germ1link$corr <- abs(germ1link$corr)

germ2link <- germ2[,c(1,2,4)]
germ2link$corr <- abs(germ2link$corr)


germAllNet$corr <- abs(germAllNet$corr)
allTypesNet$corr <- abs(allTypesNet$corr)


write.table(germ2,"afterDefence/germNet-FC2.txt",sep = "\t",row.names = F,quote = F)
write.table(germ2link,"afterDefence/germNet-FC2-linkList.txt",sep = "\t",row.names = F,quote = F,col.names = F)
