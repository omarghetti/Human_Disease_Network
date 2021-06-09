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
edges <- read.csv("Dataset/diseasome [Edges].csv", head=TRUE)
nodes <- read.csv("Dataset/diseasome [Nodes].csv", head=TRUE)
nodes <- nodes %>% select(-timeset)
edges <- edges %>% select(-timeset,-label,-id)
edges$Target <- as.integer(edges$Target)
network_graph <- graph_from_data_frame(edges, directed = TRUE,vertices=nodes)
nd3_graph <- igraph_to_networkD3(network_graph,group=nodes$X0,what = "both")
network_g <- forceNetwork(Links=nd3_graph[["links"]], 
                          Nodes=nd3_graph[["nodes"]], 
                          Source="source", 
                          Target="target",
                          NodeID="name",
                          Group="group",
                          fontSize = 8,
                          legend = TRUE,
                          zoom = TRUE,
                          radiusCalculation ="Math.sqrt(d.nodesize)+1"
                          )
v_count=vcount(network_graph)
g_count=sum(V(network_graph)$X1 == "gene")
d_count=sum(V(network_graph)$X1!="gene")
e_count=ecount(network_graph)
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
    tabPanel(
        "e_analysis",
        fluidRow(
            column(12, wellPanel(forceNetworkOutput("network_plot", width = "100%", height = "700px")))
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
    )
    tabPanel
)

# Define UI for application that draws a histogram
ui <- dashboardPage(header, sidebar, body)

# Define server logic required to draw a histogram
server <- function(input, output) {
    output$network_plot <- renderForceNetwork({
        network_g
    })

}

# Run the application 
shinyApp(ui = ui, server = server)
