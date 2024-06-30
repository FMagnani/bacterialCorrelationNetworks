# NETWORK CONNECTIVITY
library(igraph)
library(VennDiagram)

load(path_adjacencies)

#-------------------------------------------------------------------------------

# Link in comune

# Same number since they have same density
N_links <- sum(abs(adjacencies$RA_adj))

# RA-RD
both_one <- sum((adjacencies$RA_adj==1)&(adjacencies$RD_adj==1))
both_minus_one <- sum((adjacencies$RA_adj==-1)&(adjacencies$RD_adj==-1))
100*(both_one+both_minus_one)/39008

# RA-RL
both_one <- sum((adjacencies$RA_adj==1)&(adjacencies$RL_adj==1))
both_minus_one <- sum((adjacencies$RA_adj==-1)&(adjacencies$RL_adj==-1))
100*(both_one+both_minus_one)/39008

# RD-RL
both_one <- sum((adjacencies$RD_adj==1)&(adjacencies$RL_adj==1))
both_minus_one <- sum((adjacencies$RD_adj==-1)&(adjacencies$RL_adj==-1))
100*(both_one+both_minus_one)/39008

# RA-RD-RL
both_one <- sum((adjacencies$RA_adj==1)&(adjacencies$RD_adj==1)&(adjacencies$RL_adj==1))
both_minus_one <- sum((adjacencies$RA_adj==-1)&(adjacencies$RD_adj==-1)&(adjacencies$RL_adj==-1))
100*(both_one+both_minus_one)/39008

#-------------------------------------------------------------------------------

get_graph_stats <- function(adj){
  
  g <- graph_from_adjacency_matrix(
    adj,
    mode = "undirected", weighted = TRUE, diag = FALSE, 
    add.colnames = NA, add.rownames = NULL
  )
  
  comp <- components(g)
  
  print(table(table(comp$membership)))
  
  label_LCC <- names(table(comp$membership))[table(comp$membership)==max(table(comp$membership))]
  nodes_LCC <- names(comp$membership)[comp$membership==label_LCC]
  isolated_nodes <- names(table(comp$membership))[table(comp$membership)==1]

  return(list("nodes_LCC"=nodes_LCC, "isolated_nodes"=isolated_nodes))
  
}

RA_graph_stats <- get_graph_stats(adjacencies$RA_adj)
RD_graph_stats <- get_graph_stats(adjacencies$RD_adj)
RL_graph_stats <- get_graph_stats(adjacencies$RL_adj)

#-------------------------------------------------------------------------------
RA:
1    2    3    4   13 1364 
535   20    5    2    1    1 

RD:
1    2    3    4    6    8 1276 
618   21    7    1    1    1    1 

RL:
1    2    3    4    5 1345 
555   23    4    3    1    1
#-------------------------------------------------------------------------------

RA_LCC <- RA_graph_stats$nodes_LCC
RA_isolated <- RA_graph_stats$isolated_nodes
RD_LCC <- RD_graph_stats$nodes_LCC
RD_isolated <- RD_graph_stats$isolated_nodes
RL_LCC <- RL_graph_stats$nodes_LCC
RL_isolated <- RL_graph_stats$isolated_nodes

venn.diagram(
  x = list(
    RA_LCC, RD_LCC, RL_LCC
  ), 
  category.names = c("RA","RD","RL"),
  filename="LCC_venn.png",
  print.mode="percent",
  output=TRUE, imagetype="tiff", resolution=600,
  col=c("#F8766D", '#00BA38', '#619CFF'),
  fill=c("#F8766D", '#00BA38', '#619CFF')
)

venn.diagram(
  x = list(
    RA_isolated, RD_isolated, RL_isolated
  ), 
  category.names = c("RA","RD","RL"),
  filename="isolated_venn.png",
  print.mode="percent",
  output=TRUE, imagetype="tiff", resolution=600,
  col=c("#F8766D", '#00BA38', '#619CFF'),
  fill=c("#F8766D", '#00BA38', '#619CFF')
)

#-------------------------------------------------------------------------------
