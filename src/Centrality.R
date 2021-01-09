## Libraries and Graph Loading
source("Utils.R")
libraries <- c("ggraph","igraph","dplyr","readr", "DiagrammeR", "tidyverse", "Cairo", 
               "networkD3","CINNA","scales","pander")
import_libraries(libraries)
edges <- read.csv("Dataset/diseasome [Edges].csv", head=TRUE)
nodes <- read.csv("Dataset/diseasome [Nodes].csv", head=TRUE)

## Cleaning
nodes <- nodes %>% select(-timeset)
edges <- edges %>% select(-timeset,-label)

network_graph <- graph.data.frame(edges,directed = TRUE, vertices=nodes)

#Console Printing of the graph
print(network_graph, e=TRUE, v=TRUE)

#Plotting The Graph(Green nodes are genes, Blue Nodes are Diseases)
ggraph(network_graph, layout="graphopt") +
  geom_edge_fan(colour="grey15",show.legend = TRUE) +
  geom_node_point(fill= ifelse(nodes$X0 == "gene", "#00FF00", "#0D98BA"), 
                  shape=23, col="grey15", show.legend = TRUE) +
  scale_size_continuous(range=c(1, 10)) +
  theme_graph(base_size = 11, base_family = "sans") +
  ggtitle("Human Disease Network") 

#Exploratory Analysis
print(paste("Vertices Number:", vcount(network_graph)))
print(paste("Genes Number:", sum(V(network_graph)$X1 == "gene")))
print(paste("Diseases Number:", sum(V(network_graph)$X1!="gene")))
print(paste("Edges Number:", ecount(network_graph)))

##Centrality
#Degree Centrality
network_degree <- centr_degree(network_graph,mode="all",normalized = TRUE)

  

