import_libraries <- function(libraries_list) {
  if(!require("easypackages")) {
    install.packages("easypackages")
  }
  library(easypackages)
  packages(libraries_list)
  libraries(libraries_list)
}

Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

histogram_plot <- function(mapping, scale, x_lab, y_lab, title){
  ggplot(mapping=aes(x=mapping)) + 
    geom_histogram(col="black", fill="#8B0000") + 
    stat_bin(aes(y=..count.., label=..count..), geom="text", vjust=-.5, size=3) +
    scale_y_log10(labels = trans_format("log10", math_format(expr = 10^.x))) +
    scale_x_continuous(x_lab, scale) +
    labs(x=x_lab, y=y_lab) +
    ggtitle(title) +
    theme_bw()
}

plot_graph <- function(graph, layout, measure, node_s1, node_s2, node_s3, label, label_s1, title){
  ggraph(graph, layout=layout) + #  "graphopt""
    geom_edge_fan(colour = "gray66") +
    geom_node_point(aes(fill="yellow", 
                        size=ifelse(measure > node_s1, 30,
                                    ifelse(measure > node_s2,  8, 
                                           ifelse(measure > node_s3,  3, 1)))),
                    shape=21, col="grey25", show.legend = FALSE) +
    geom_node_text(aes(label = ifelse(measure > label_s1, 
                                      as.character(label), NA), size=2), 
                   show.legend = FALSE) +
    scale_size_continuous(range=c(1, 10)) +
    theme_graph(base_size = 11, base_family = "serif") +
    ggtitle(title) 
}

plot_clustering_graph_with_legend <- function(graph, layout, fill, edge_width, title){
  ggraph(graph, layout=layout) +
    geom_edge_fan(aes(width=edge_width), colour = "gray66", show.legend = FALSE) +
    geom_node_point(aes(fill=fill), shape=21, size=3) +
    scale_fill_manual(values = rainbow(length(table(fill)))) +
    scale_edge_width_continuous(range=c(0.2,0.9)) +
    labs(fill="Clusters") +
    guides(fill = guide_legend(title = "Cluster Colors", title.position = "top"), col = guide_legend(ncol = 5)) +
    theme_graph(base_size = 11, base_family = "serif") +
    theme(legend.position = "bottom", legend.text=element_text(size=8)) +
    ggtitle(title) 
}

confusion_matrix = function(list1, list2){
  finalList = union(list1,list2)
  list1[list1=="Multiple"] = list2[list1=="Multiple"]
  list1[list1=="Unclassified"] = list2[list1=="Unclassified"]
  new_table = table(factor(list1,finalList),factor(list2,finalList))
  return(confusionMatrix(new_table))
}

