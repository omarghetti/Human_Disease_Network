## Libraries and Graph Loading
source("Shared_Functions.R")
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

#Centrality
#Degree Centrality
network_degree = centr_degree(network_graph,mode="all",normalized = TRUE)

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

paste("Nodes: ", vcount(network_graph_no_genes))
paste("Edges: ", ecount(network_graph_no_genes))

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
clustering_result <- nodes %>% filter(X1 != "gene") %>% select(-X0)
#Girvan-Newmann
girvan_newmann = cluster_edge_betweenness(network_graph_no_genes,directed=FALSE,
                                          weights=E(network_graph_no_genes)$weight) 
print(paste("GN Number Of Communities:",max(girvan_newmann$membership)))
girvannewmann_df = label_cluster(girvan_newmann)
clustering_result$girvannewmann <- NA
for (n in 1:length(girvan_newmann$membership))
{
  clustering_result$girvannewmann[n] = girvannewmann_df$name[girvan_newmann$membership[n]]
}
plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", nodes_results$girvannewmann,
                         E(network_graph_no_genes)$weight, "Girvan Newman")

##Louvain
louvain_clustering = cluster_louvain(network_graph_no_genes, weights = E(network_graph_no_genes)$weights)

print(paste("Louvain Number Of Communities:",max(louvain_clustering$membership)))
louvain_df = label_cluster(louvain_clustering)
clustering_result$Louvain <- NA
for (n in 1:length(louvain_clustering$membership))
{
  clustering_result$Louvain[n] = louvain_df$name[louvain_clustering$membership[n]]
}
plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", nodes_results$Louvain,
                                  E(network_graph_no_genes)$weight, "Louvain")
##Fastgreedy
fastgreedy_clusters  = cluster_fast_greedy(network_graph_no_genes, modularity = TRUE,
                                           weights = E(network_graph_no_genes)$weight)

print(paste("Fastgreedy communities: ", max(fastgreedy_clusters$membership)))
fgreedy_df <- label_cluster(fastgreedy_clusters)
clustering_result$Fastgreedy <- NA
for (n in 1:length(fastgreedy_clusters$membership))
{
  clustering_result$Fastgreedy[n] <- fgreedy_df$name[fastgreedy_clusters$membership[n]]
}
plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", nodes_results$Fastgreedy,
                         E(network_graph_no_genes)$weight, "Fastgreedy")
##Leiden
mat = as_adj(network_graph_no_genes, type = "both", attr = "weight")
leiden_clustering <- leiden(mat, resolution_parameter = 86)
print(max(leiden_clustering))
leiden_clustering <- make_clusters(network_graph_no_genes, membership = leiden_clustering)
leiden_df <- label_cluster(leiden_clustering)
clustering_result$Leiden = NA
for (n in 1:length(leiden_clustering$membership))
{
  clustering_result$Leiden[n] <- leiden_df$name[leiden_clustering$membership[n]]
}
plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", clustering_result$Leiden,
                         E(network_graph_no_genes)$weight, "Leiden")
##Markov Clustering
mat = as_adj(network_graph_no_genes,type = "both",attr = "weight")
markov_clustering = mcl(mat,addLoops = FALSE, inflation=4, allow1 = TRUE)
print(paste("Markov Communities:",markov_clustering$K))
markov_df <- as.data.frame(matrix(1:(markov_clustering$K), nrow = markov_clustering$K, 
                                          dimnames = list(NULL, "id")))
markov_df$name <- NA
markov_df$id <- unique(markov_clustering$Cluster)
for (id in markov_df$id)
{
  labels <- V(network_graph_no_genes)[markov_clustering$Cluster == id]$X1
  labels <- labels[labels != "Multiple"]
  
  value <- which.max(unlist(table(labels)))
  
  if (length(value) != 0)
  {
    markov_df$name[markov_df$id == id] <- names(value)[1]
  }
  else
  {
    markov_df$name[markov_df$id == id] <- "Multiple"
  }
}
clustering_result$Markov <- NA
for (n in 1:length(markov_clustering$Cluster))
{
  clustering_result$Markov[n] <- markov_df$name[markov_clustering$Cluster[n] == markov_df$id]
}
plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", nodes_results$Markov,
                         E(network_graph_no_genes)$weight, "Markov Clusters")
##Leading EigenVector
lead_eigen <- cluster_leading_eigen(network_graph_no_genes, weights = E(network_graph_no_genes)$weight)
print(paste("Leading Eigen Communities: ", max(lead_eigen$membership)))
lead_eigen_df <- label_cluster(lead_eigen)
clustering_result$Lead_eigen = NA
for (n in 1:length(lead_eigen$membership))
{
  clustering_result$Lead_eigen[n] <- lead_eigen_df$name[lead_eigen$membership[n]]
}
plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", nodes_results$Lead_eigen,
                         E(network_graph_no_genes)$weight, "Leading eigenvector Clusters")
##Ground-Truth for comparison
groundtruth_clusters = unique(V(network_graph_no_genes)$X1)
i = 1
for (value in groundtruth_clusters){
  V(network_graph_no_genes)[V(network_graph_no_genes)$X1 == value]$color = i
  i = i + 1
}

plot_clustering_graph_with_legend(network_graph_no_genes, "graphopt", nodes$X1[nodes$X0 == "disease"], 
                         E(network_graph_no_genes)$weight, "Groundtruth Clustering")

#Comparing Algorithm Results
#Purity
install.packages("e1071")
p_frame = data.frame(matrix(ncol = 23, nrow = 7))
colnames(p_frame) = c(groundtruth_clusters,"Purity Measure")
rownames(p_frame) = c("newmann","louvain","fgreedy","leiden","markov","lead_eigen","avg")

for(i in 4:length(clustering_result)){
  cmatrix = confusion_matrix(clustering_result$X1,clustering_result[,i])
  p_frame[i-3, ]=unname(c(diag(cmatrix$table/rowSums(cmatrix$table)),
                          cmatrix$overall[1])) 
}
p_frame$Multiple = NULL
p_frame$Unclassified = NULL
p_frame[7, ] = colMeans(p_frame[c(1:6), ])
p_frame = as.data.frame(t(p_frame))
panderOptions('table.split.table',8*15)
pander(format(p_frame))
#F-Measure
f_frame = data.frame(matrix(ncol = 1, nrow = 6))
for(i in 4:length(clustering_result)){
  cmatrix = confusion_matrix(clustering_result$X1,clustering_result[,i])
  f_frame[i-3, ] = sum(cmatrix$byClass[ ,7], na.rm = TRUE) / 22
}
colnames(f_frame) = c("F-Measure")
rownames(f_frame) = c("newmann","louvain","fgreedy","leiden","markov","lead_eigen")
pander(f_frame)

#Normalized Mutual Information
mutualinfo_df <- data.frame(matrix(ncol = 1, nrow = 6))
colnames(mutualinfo_df) <- c("NMI")
rownames(mutualinfo_df) <- c("newmann", "louvain", "fgreedy", "leiden",
                      "markov", "lead_eigen")

for (i in 1:nrow(mutualinfo_df))
{
  groundtruth = clustering_result[, 3]
  clusters = clustering_result[, 3 + i]
  groundtruth[groundtruth == "Multiple"] <- clusters[groundtruth == "Multiple"]
  groundtruth[groundtruth == "Unclassified"] <- clusters[groundtruth == "Unclassified"]
  mutualinfo_df[i, ] <- NMI(groundtruth, clusters)
}
pander(format(mutualinfo_df))

#Adjusted Rand Index
install.packages("aricode")
rand_df = data.frame(matrix(ncol=1,nrow=6))
colnames(rand_df) = c("Adjusted Rand Index")
rownames(rand_df) = c("newmann", "louvain", "fgreedy", "leiden",
                      "markov", "lead_eigen")
for (i in 1:nrow(rand_df))
{
  groundtruth <- clustering_result[, 3]
  clusters <- clustering_result[, 3 + i]
  groundtruth[groundtruth == "Multiple"] <- clusters[groundtruth == "Multiple"]
  groundtruth[groundtruth == "Unclassified"] <- clusters[groundtruth == "Unclassified"]
  rand_df[i, ] <- ARI(groundtruth, clusters)
}
pander(format(rand_df))
