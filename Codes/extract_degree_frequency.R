rm(list = ls())

library('stats4')
library('igraph')
library('Matrix')
library('lattice')
#library('sna')  ## do not use 'sna' package and 'igraph' packages simultaneous , because it si not compatible.
library('MASS')
library('base')
library('ggplot2')



#data<-read.table('WebGoogleClean.txt')
#edgeset=cbind((data$V1+1),(data$V2+1))


#data<-read.table('Cit_HepPhClean.txt')
#edgeset=cbind(data$V1,data$V2)


#data<-read.table('AstroPhClean.txt')
#edgeset=cbind(data$V1,data$V2)


#############  for indegree ego_gplus network  ######################################

data<-read.table('bio-SC-HT_edges.txt',numerals = c("no.loss"))
all_vertices_array=append(data$V1,data$V2)
#all_vertices_array=data$V2
table_array=table(all_vertices_array)
new_table_array=as.data.frame(table_array)
frequency_table=table(new_table_array$Freq)
final_frequency_table=as.data.frame(frequency_table)

write.csv(final_frequency_table,'final_degree_frequency_bio-SC-HT.csv',sep=',',row.names = FALSE)

#######################################################################################


    # data<-read.table('gplus_combined.txt',numerals = c("no.loss"))
    # #edgeset=cbind((data$node_1+1),(data$node_2+1))
    # edgeset=cbind((data$V1),(data$V2))
    # 
    # ###  for directed graph
    # g1=graph_from_edgelist(edgeset,dir=TRUE)
    # 
    # ###  for undirected graph
    # #g1=graph_from_edgelist(edgeset,dir=FALSE)
    # 
    # 
    # 
    # g1<-simplify(g1,remove.multiple = TRUE, remove.loops = TRUE)
    # 
    # 
    # ###  for directed graph
    # node_degree=degree(g1,v=V(g1),mode='in')  
    # 
    # ###  for undirected graph
    # #node_degree=degree(g1)
    # 
    # 
    # table_freq=table(node_degree)
    # 
    # degree_frequency_table=as.data.frame(table_freq)
    # 
    # 
    # write.csv(degree_frequency_table,'final_indegree_frequency_gplus.csv',sep=',',row.names = FALSE)
    # 
