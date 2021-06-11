#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
library(shinydashboard)
library(shiny)
source("Shared_Functions.R")
libraries <- c("ggraph","igraph","ggiraph","dplyr","readr", "DiagrammeR", "tidyverse", "Cairo", 
               "networkD3","CINNA","scales","pander", "plotly")
import_libraries(libraries)
libraries = c("scales","gridExtra","leiden","pander","MCL","caret","aricode")
import_libraries(libraries)
edges <- read.csv("Dataset/diseasome [Edges].csv", head=TRUE)
nodes <- read.csv("Dataset/diseasome [Nodes].csv", head=TRUE)
nodes <- nodes %>% select(-timeset)
edges <- edges %>% select(-timeset,-label,-id)
edges$Target <- as.integer(edges$Target)
network_graph <- graph_from_data_frame(edges, directed = TRUE,vertices=nodes)
nd3_graph <- igraph_to_networkD3(network_graph,group=nodes$X1,what = "both")
network_g <- forceNetwork(Links=nd3_graph[["links"]], 
                          Nodes=nd3_graph[["nodes"]], 
                          Source="source", 
                          Target="target",
                          NodeID="name",
                          Group="group",
                          fontSize = 12,
                          zoom = TRUE,
                          legend = TRUE,
                          opacity=0.8,
                          linkDistance = 1,
                          width = 700,
                          height = 1000,
                          radiusCalculation ="Math.sqrt(d.nodesize)+1"
                          )
v_count=vcount(network_graph)
g_count=sum(V(network_graph)$X1 == "gene")
d_count=sum(V(network_graph)$X1!="gene")
e_count=ecount(network_graph)
network_degree = centr_degree(network_graph,mode="all",normalized = TRUE)
b_cetr <- betweenness(network_graph,v=V(network_graph),directed=TRUE,normalized=TRUE)
c_centr <- closeness(network_graph_no_genes, vids=V(network_graph_no_genes),mode="all",normalized=TRUE)
eigenv <- eigen_centrality(network_graph_no_genes,directed=FALSE)
prank <- page_rank(network_graph_no_genes,vids=V(network_graph_no_genes),directed = FALSE)
for(i in 1:length(edges[,3])){
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

new_graph = graph_from_data_frame(edges,directed=FALSE,vertices=nodes)
network_graph_no_genes = induced_subgraph(new_graph,which(nodes$X1 != "gene"))
network_graph_no_genes = igraph::simplify(network_graph_no_genes)
nodes_no_genes = nodes[nodes$X1!="gene", ]
nd3_graph_new = igraph_to_networkD3(network_graph_no_genes,group=nodes_no_genes$X1,what = "both")
network_g_new <- forceNetwork(Links=nd3_graph_new[["links"]], 
                          Nodes=nd3_graph_new[["nodes"]], 
                          Source="source", 
                          Target="target",
                          Value = "value",
                          NodeID="name",
                          Group="group",
                          fontSize = 12,
                          zoom = TRUE,
                          legend = TRUE,
                          opacity=0.8,
                          linkDistance = 1,
                          width = 700,
                          height = 1000,
                          radiusCalculation ="Math.sqrt(d.nodesize)+1"
)
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
clustering_result <- nodes %>% filter(X1 != "gene") %>% select(-X0)
girvan_newmann = cluster_edge_betweenness(network_graph_no_genes,directed=FALSE,
                                          weights=E(network_graph_no_genes)$weight) 
print(paste("GN Number Of Communities:",max(girvan_newmann$membership)))
girvannewmann_df = label_cluster(girvan_newmann)
clustering_result$girvannewmann <- NA
for (n in 1:length(girvan_newmann$membership))
{
    clustering_result$girvannewmann[n] = girvannewmann_df$name[girvan_newmann$membership[n]]
}
louvain_clustering = cluster_louvain(network_graph_no_genes, weights = E(network_graph_no_genes)$weights)

print(paste("Louvain Number Of Communities:",max(louvain_clustering$membership)))
louvain_df = label_cluster(louvain_clustering)
clustering_result$Louvain <- NA
for (n in 1:length(louvain_clustering$membership))
{
    clustering_result$Louvain[n] = louvain_df$name[louvain_clustering$membership[n]]
}
fastgreedy_clusters  = cluster_fast_greedy(network_graph_no_genes, modularity = TRUE,
                                           weights = E(network_graph_no_genes)$weight)

print(paste("Fastgreedy communities: ", max(fastgreedy_clusters$membership)))
fgreedy_df <- label_cluster(fastgreedy_clusters)
clustering_result$Fastgreedy <- NA
for (n in 1:length(fastgreedy_clusters$membership))
{
    clustering_result$Fastgreedy[n] <- fgreedy_df$name[fastgreedy_clusters$membership[n]]
}
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
lead_eigen <- cluster_leading_eigen(network_graph_no_genes, weights = E(network_graph_no_genes)$weight)
print(paste("Leading Eigen Communities: ", max(lead_eigen$membership)))
lead_eigen_df <- label_cluster(lead_eigen)
clustering_result$Lead_eigen = NA
for (n in 1:length(lead_eigen$membership))
{
    clustering_result$Lead_eigen[n] <- lead_eigen_df$name[lead_eigen$membership[n]]
}
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

#Define Application Header
header <- dashboardHeader(title = HTML("HUMAN DISEASE NETWORK ANALYSIS"))
sidebar <- dashboardSidebar(
    width = 200,
    sidebarMenu(
        id='sidebar',
        style="position: relative, overflow: visible",
        menuItem(
            "Exploratory Analysis",
            tabName = "e_analysis",
            icon = icon('dashboard'),
            badgeColor = 'green'
        ),
        menuItem(
            "Centrality Analysis",
            tabName = "ce_analysis",
            icon = icon('dashboard'),
            badgeColor = 'green'
        ),
        menuItem(
            "Clustering Analysis",
            tabName = "cl_analysis",
            icon = icon('dashboard'),
            badgeColor = 'green'
        )
    )
)
body <- dashboardBody(
    tabItems(
    tabItem(
        "e_analysis",
        fluidRow(
            column(12, wellPanel(forceNetworkOutput("network_plot", width = "100%", height = "1000px")))
        ),
        h1(paste0("Exploratory Analysis")),
        fluidRow(
          valueBox(v_count,"Vertices Number",color = "aqua"),
          valueBox(g_count,"Genes Number", color = "green")
        ),
        fluidRow(
          valueBox(d_count, "Diseases Number", color = "red"),
          valueBox(e_count, "Edges Number", color = "orange")
        )
    ),
    tabItem(
        "ce_analysis",
        h1(paste0("Centrality Analysis")),
        fluidRow(
            column(12,
                   wellPanel(plotlyOutput("Degreecentr"))),
            column(12,
                   wellPanel(plotlyOutput("Bcentr"))),
            column(12,
                   wellPanel(plotlyOutput("Ccentr"))),
            column(12,
                   wellPanel(plotlyOutput("Pcentr"))),
            column(12,
                   wellPanel(plotlyOutput("Ecentr"))),
        )),
    tabItem("cl_analysis",
            h1(paste0("Clustering Analysis")),
            h2(paste0("Network Prepared for Clustering Analysis")),
            fluidRow(
                column(12, wellPanel(forceNetworkOutput("new_network_plot", width = "100%", height = "1000px")))
            ),
            h1(paste0("Purity Comparison")),
            fluidRow(
                column(12, tableOutput('clustering_table'))
            ))
    )
)

# Define UI for application that draws a histogram
ui <- dashboardPage(header, sidebar, body)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$network_plot <- renderForceNetwork({
        network_g
    })
    
    output$new_network_plot <- renderForceNetwork({
        network_g_new
    })
    
    output$Degreecentr <- renderPlotly({
        ggplotly(histogram_plot(network_degree$res,seq(0,180,by=20),"n_degree","freq","Degree Centrality"))
    })
    
    output$Bcentr <- renderPlotly({
        ggplotly(histogram_plot(b_cetr, c(0,1),"betweenness","freq","Betweenness Centrality"))
    })
    
    output$Ccentr <- renderPlotly({
        histogram_plot(c_centr,c(0,1), "Closeness","Freq","Closeness")
    })
    
    output$Pcentr <- renderPlotly({
        histogram_plot(eigenv$vector, c(0,1),"eigen_N","Freq","Eigenvector")
    })
    
    output$Ecentr <- renderPlotly({
        histogram_plot(prank$vector,c(0,1),"Pagerank","Freq","Pagerank")
    })
    
    output$clustering_table <-renderTable({
        p_frame
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
