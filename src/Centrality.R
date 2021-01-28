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
#plot_graph(network_graph, "dh", b_cetr, 0.14, 0.01, 0.0001, nodes$label, 0.14, 
#          "Human Disease Network Betweenness Centrality")

#Closeness
c_centr <- closeness(network_graph, vids=V(network_graph),mode="all",normalized=TRUE)
histogram_plot(c_centr,c(0,1), "Closeness","Freq","Closeness")
#plot_graph(network_graph, "dh", c_centr, 0.14, 0.01, 0.0001, nodes$label, 0.14, 
#          "Human Disease Network Betweenness Centrality")
#EigenVector
eigenv <- eigen_centrality(network_graph,directed=TRUE)
histogram_plot(eigenv$vector, c(0,1),"eigen_N","Freq","Eigenvector")

#PageRank
prank <- page_rank(network_graph,vids=V(network_graph),directed = TRUE)
histogram_plot(prank$vector,c(0,1),"Pagerank","Freq","Pagerank")

#Clustering Coefficient
g_trans <- transitivity(network_graph, type="global",isolates="zero")
l_trans <- transitivity(network_graph, type="local",isolates="zero")
transitivity_recap <- data.frame(matrix(ncol=1, nrow=2))
colnames(transitivity_recap) <- "Values"
rownames(transitivity_recap) <- c("Local Transitivity Average","Global Transitivity")
transitivity_recap[1,1] <- mean(l_trans)
transitivity_recap[2,1] <- g_trans
pander(transitivity_recap,digits=4,justify="center")

#Moving To a Non-Directed Graph to facilitate the Clustering Analysis, Recalculating Weights
for(i in 1:length(edges[, 5])){
  e = edges[i, ]
  source_node = V(network_graph)[as.character(e$Source)]
  target_node = V(network_graph)[as.character(e$Target)]
  neigh_source_node = neighborhood(network_graph,nodes=source_node)[[1]]
  neigh_target_node = neighborhood(network_graph,nodes=target_node)[[1]]
  neigh_source_node = as.numeric(neigh_source_node[neigh_source_node&X1 == "gene"])
  neigh_target_node = as.numeric(neigh_target_node[neigh_target_node&X1 == "gene"])
  new_weight = length(intersect(neigh_source_node,neigh_target_node))
  edges[i, ]$weight = new_weight
}

new_graph = graph.data.frame(edges,directed=FALSE,vertices=nodes)

#REMOVING GENES
network_graph_no_genes = induced_subgraph(new_graph,which(nodes$X1 != "gene"))
network_graph_no_genes = igraph::simplify(network_graph_no_genes)
ggraph(network_graph_no_genes, layout="graphopt") +
  geom_edge_fan(aes(width=E(network_graph_no_genes)$weight), colour= "gray66",show.legend = FALSE) +
  geom_node_point(fill= "#00FF00",
                  shape=21, col="gray15", show.legend = FALSE) +
  scale_edge_width_continuous(range=c(0.2,0.9)) +
  scale_size_continuous(range=c(1, 10)) +
  theme_graph(base_size = 11, base_family = "sans") +
  ggtitle("Human Disease Network With Genes Removed") 

paste("N° di nodi del grafo finale: ", vcount(network_graph_no_genes))
paste("N° di archi del grafo finale: ", ecount(network_graph_no_genes))

#New Graph Network Analysis
#Degree Centrality
network_degree <- centr_degree(network_graph_no_genes,mode="all",normalized = TRUE)

histogram_plot(network_degree$res,seq(0,180,by=20),"n_degree","freq","Degree Centrality")

#Top Nodes by Degree
V(network_graph)$label[order(network_degree$res, decreasing = TRUE)][1:10]

#Bottom Nodes by Degree
V(network_graph)$label[order(network_degree$res, decreasing = FALSE)][1:10]

#Betweenness Centrality
b_cetr <- betweenness(network_graph_no_genes,v=V(network_graph_no_genes),directed=FALSE,normalized=TRUE)
histogram_plot(b_cetr, c(0,1),"betweenness","freq","Betweenness Centrality") 

#Closeness
c_centr <- closeness(network_graph_no_genes, vids=V(network_graph_no_genes),mode="all",normalized=TRUE)
histogram_plot(c_centr,c(0,1), "Closeness","Freq","Closeness")
#plot_graph(network_graph, "dh", c_centr, 0.14, 0.01, 0.0001, nodes$label, 0.14, 
#          "Human Disease Network Betweenness Centrality")
#EigenVector
eigenv <- eigen_centrality(network_graph_no_genes,directed=FALSE)
histogram_plot(eigenv$vector, c(0,1),"eigen_N","Freq","Eigenvector")

#PageRank
prank <- page_rank(network_graph_no_genes,vids=V(network_graph_no_genes),directed = FALSE)
histogram_plot(prank$vector,c(0,1),"Pagerank","Freq","Pagerank")

#Clustering Coefficient
g_trans <- transitivity(network_graph_no_genes, type="global",isolates="zero")
l_trans <- transitivity(network_graph_no_genes, type="local",isolates="zero")
transitivity_recap <- data.frame(matrix(ncol=1, nrow=2))
colnames(transitivity_recap) <- "Values"
rownames(transitivity_recap) <- c("Local Transitivity Average","Global Transitivity")
transitivity_recap[1,1] <- mean(l_trans)
transitivity_recap[2,1] <- g_trans
pander(transitivity_recap,digits=4,justify="center")

##Clustering Section
libraries = c("scales","gridExtra","leiden","pander","MCL","caret","aricode")
import_libraries(libraries)

##Clustering Algorithms
#Function to label clusters
label_cluster <- function(clusters)
{
  cluster_df <-
    as.data.frame(matrix(
      1:length(clusters),
      nrow = length(clusters),
      dimnames = list(NULL, "id")
    ))
  
  for (c in 1:length(clusters))
  {
    labels <- V(network_graph_no_genes)[clusters[[c]]]$X1
    labels <- labels[which(labels != "Multiple" & labels != "Unclassified")]
    
    if (length(labels) != 0)
    {
      value <- which.max(unlist(table(labels)))
      cluster_df$name[cluster_df$id == c] <- names(value)[1]
    } else
    {
      cluster_df$name[cluster_df$id == c] <- V(network_graph_no_genes)[clusters[[c]]]$X1[1]
    }
  }
  
  return(cluster_df)
}

#defining data-frame to compare clusters
nodes_results <- nodes %>% filter(X1 != "gene") %>% select(-X0)
nodes_results$girvannewmann <- NA
nodes_results$Louvain <- NA
nodes_results$Fastgreedy <- NA
nodes_results$Markov <- NA
nodes_results$Leiden <- NA
nodes_results$Label_prop <- NA
nodes_results$Lead_eigen <- NA

#Girvan-Newmann
girvan_newmann = cluster_edge_betweenness(network_graph_no_genes,directed=FALSE,
                                          weights=E(network_graph_no_genes)$weight) 
print(paste("GN Number Of Communities:",max(girvan_newmann$membership)))
girvannewmann_df = label_cluster(girvan_newmann)
for (n in 1:length(girvan_newmann$membership))
{
  nodes_results$girvannewmann[n] = girvannewmann_df$name[girvan_newmann$membership[n]]
}
plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", nodes_results$girvannewmann,
                         E(network_graph_no_genes)$weight, "Girvan Newman")

##Louvain
louvain_clustering = cluster_louvain(network_graph_no_genes, weights = E(network_graph_no_genes)$weights)

print(paste("Louvain Number Of Communities:",max(louvain_clustering$membership)))
louvain_df = label_cluster(louvain_clustering)
for (n in 1:length(louvain_clustering$membership))
{
  nodes_results$Louvain[n] = louvain_df$name[louvain_clustering$membership[n]]
}
plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", nodes_results$Louvain,
                                  E(network_graph_no_genes)$weight, "Louvain")
##Fastgreedy
fastgreedy_clusters  = cluster_fast_greedy(network_graph_no_genes, modularity = TRUE,
                                           weights = E(network_graph_no_genes)$weight)

print(paste("N° di communities: ", max(fastgreedy_clusters$membership)))
fgreedy_df <- label_cluster(fastgreedy_clusters)
for (n in 1:length(fastgreedy_clusters$membership))
{
  nodes_results$Fastgreedy[n] <- fgreedy_df$name[fastgreedy_clusters$membership[n]]
}