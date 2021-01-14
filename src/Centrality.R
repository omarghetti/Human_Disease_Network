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

histogram_plot(network_degree$res,seq(0,180,by=20),"n_degree","freq","Degree Centrality")

#Directed Graph: In and Out Degree
odegree <- centr_degree(network_graph, mode="out",normalized=TRUE)
histogram_plot(odegree$res, seq(0,20,by=20),"outdegree","freq","OutDegree")
idegree <- centr_degree(network_graph, mode="in",normalized=TRUE)
histogram_plot(idegree$res, seq(0,20,by=20),"indegree","freq","InDegree")

#Top Nodes by Degree
V(network_graph)$label[order(network_degree$res, decreasing = TRUE)][1:10]

#Bottom Nodes by Degree
V(network_graph)$label[order(network_degree$res, decreasing = FALSE)][1:10]

#Betweenness Centrality
b_cetr <- betweenness(network_graph,v=V(network_graph),directed=TRUE,normalized=TRUE)
histogram_plot(b_cetr, c(0,1),"betweenness","freq","Betweenness Centrality") 
plot_graph(network_graph, "dh", b_cetr, 0.14, 0.01, 0.0001, nodes$label, 0.14, 
           "Human Disease Network Betweenness Centrality")

#Closeness
c_centr <- closeness(network_graph, vids=V(network_graph),mode="all",normalized=TRUE)