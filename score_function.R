## we are using R graph
#install.packages("igraph")
library(igraph)



DElist1 <- read.delim("2data-cb-easyDElist.txt",sep = "\t",header = T)
colnames(DElist1)[1] <- "Name"
DElist1 <- DElist1[,1]

cellnet <- read.delim("Mouse_Big_GRN.txt",sep = " ",header = T)
cellnet <- cellnet[,c(2,1,3:6)]

cellnetGraph<-graph_from_data_frame(cellnet,directed=T)

DElist1 = as.character(DElist1)

# g : igraph object, sub-network
# geneList : a vector of biological process genes
# iter: number of iteration for z-score
Panir_subgraph <- function(g,geneList , iter = 1000 ) {
  
  if( is.igraph(g) == FALSE )
    stop("g should be igraph object")
  if( is.character(geneList) )
    stop("gene list shold be character")
  
  subGraph_geneNames <- intersect(geneList,V(g)$name)
  
  if( length(subGraph_geneNames) == 0 )
    stop(" there is no intersect between vertex of g and geneList")
  
  subg = induced_subgraph(g,subGraph_geneNames)
  
  intra_interaction = sum(degree(subg,mode="total" , loops = FALSE)) / 2
  
  inter_interaction = sum(degree(g, v = subGraph_geneNames , mode="total" , loops = FALSE)) - 2*intra_interaction
  
  m = intra_interaction / ( inter_interaction + intra_interaction)
  
  n = length(geneList)
  vertexNames = V(g)
  
  random_scores = lapply(1:iter, function(x){ 
    tmp_names = sample(vertexNames , n )
    tmp_sub = induced_subgraph(g,tmp_names)
    
    intra_interaction_tmp = sum(degree(tmp_sub,mode="total" , loops = FALSE)) / 2
    
    inter_interaction_tmp = sum(degree(g, v = tmp_names , mode="total" , loops = FALSE)) - 2*intra_interaction_tmp
    
    m_tmp = intra_interaction_tmp / ( inter_interaction_tmp + intra_interaction_tmp)
    return(m_tmp)
    })
  
  random_scores = unlist(random_scores)
  
  mean_scores = mean(random_scores)
  sd_scores = sd(random_scores)
  
  z =  ( m- mean_scores) / sd_scores
  print( 2*pnorm(-abs(z)))
}
