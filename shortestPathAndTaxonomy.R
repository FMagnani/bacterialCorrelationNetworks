library(ggplot2)
library(data.table)
library(gridExtra)
library(ggpubr)
library(igraph)

taxtab <- as.data.frame(fread(path_taxonomy_table.csv))

load(path_adjacencies)
species_csv <- as.data.frame(fread(path_species))

RA_adj <- adjacencies$RA_adj
RD_adj <- adjacencies$RD_adj
RL_adj <- adjacencies$RL_adj

# Take only positive links
RA_adj[RA_adj<0] <- 0
RD_adj[RD_adj<0] <- 0
RL_adj[RL_adj<0] <- 0

# Graphs
g_RA <- graph_from_adjacency_matrix(
  as.matrix(RA_adj), 
  mode = "undirected", weighted = TRUE, diag = FALSE, 
  add.colnames = NA, add.rownames = NULL
)
g_RD <- graph_from_adjacency_matrix(
  as.matrix(RD_adj), 
  mode = "undirected", weighted = TRUE, diag = FALSE, 
  add.colnames = NA, add.rownames = NULL
)
g_RL <- graph_from_adjacency_matrix(
  as.matrix(RL_adj), 
  mode = "undirected", weighted = TRUE, diag = FALSE, 
  add.colnames = NA, add.rownames = NULL
)

graph_to_lcc <- function(g){
  
  comps <- components(g)
  lcc_index <- which(comps$csize==max(comps$csize))  
  if(length(lcc_index)!=1){
    lcc_index <- lcc_index[1]
  }
  nodes <- V(g)$name
  lcc_nodes <- nodes[comps$membership==lcc_index]
  
  lcc_g <- induced_subgraph(graph=g, lcc_nodes)
  
  return(lcc_g)
  
}

# LCC
RA_lcc <- graph_to_lcc(g_RA)
RD_lcc <- graph_to_lcc(g_RD)
RL_lcc <- graph_to_lcc(g_RL)

# Shortest Paths in LCC
RA_lcc_shortestPaths <- distances(RA_lcc, weights=NA)
RD_lcc_shortestPaths <- distances(RD_lcc, weights=NA)
RL_lcc_shortestPaths <- distances(RL_lcc, weights=NA)

RA_mean <- mean(RA_lcc_shortestPaths)
RD_mean <- mean(RD_lcc_shortestPaths)
RL_mean <- mean(RL_lcc_shortestPaths)

# Now rename the distance matrix
rownames_of_matrix <- function(mat){
  speciesIds <- c()
  for(n in rownames(mat)){
    id <- species_csv[species_csv$names==n, "ncbi_ids"]
    speciesIds <- c(speciesIds, id)
  }
  return(speciesIds)
}

ids_RA <- rownames_of_matrix(RA_lcc_shortestPaths)
ids_RD <- rownames_of_matrix(RD_lcc_shortestPaths)
ids_RL <- rownames_of_matrix(RL_lcc_shortestPaths)

rownames(RA_lcc_shortestPaths) <- ids_RA
rownames(RD_lcc_shortestPaths) <- ids_RD
rownames(RL_lcc_shortestPaths) <- ids_RL
colnames(RA_lcc_shortestPaths) <- ids_RA
colnames(RD_lcc_shortestPaths) <- ids_RD
colnames(RL_lcc_shortestPaths) <- ids_RL

#-------------------------------------------------------------------------------

phylum_dict <- data.frame(
  "ids" = c(1224, 1239, 976, 201174),
  "names" = c('Pseudomonadota','Bacillota','Bacteroidota','Actinomycetota')
)

class_dict <- data.frame(
  "ids" = c(
    28211,1236,28216,117743,203490,186801,3031449,909932,91061,1760,526524,84998,3031852,200643
  ),
  "names" = c(
    'Alphaproteobacteria', 'Gammaproteobacteria', 'Betaproteobacteria', 'Flavobacteriia', 'Fusobacteriia', 'Clostridia', 'Desulfovibrionia', 'Negativicutes', 'Bacilli', 'Actinomycetes', 'Erysipelotrichia', 'Coriobacteriia', 'Epsilonproteobacteria', 'Bacteroidia'
  )
)

order_dict <- data.frame(
  "ids" = c(186802,171549,80840,186826,91347,72274,2887326,204455,200644,204457,135614),
  "names" = c('Eubacteriales','Bacteroidales','Burkholderiales','Lactobacillales','Enterobacterales','Pseudomonadales','Moraxellales','Rhodobacterales','Flavobacteriales','Sphingomonadales','Xanthomonadales')
)

#-------------------------------------------------------------------------------

get_distributions <- function(mat, tax_lvl){
  
  group_names <- unique(taxtab[, tax_lvl])
  
  results_inside <- c()
  results_cross <- c()
  
  groups_inside <- c()
  groups_cross <- c()
  
  for(g in group_names){
    
    inside_species <- taxtab[taxtab[,tax_lvl]==g, "species"]
    outside_species <- taxtab[taxtab[,tax_lvl]!=g, "species"]
    
    inside_dists <- mat[
      rownames(mat) %in% inside_species, 
      colnames(mat) %in% inside_species
    ]
    
    cross_dists <- mat[
      rownames(mat) %in% inside_species, 
      colnames(mat) %in% outside_species
    ]
    
    results_inside <- c(results_inside, as.vector(inside_dists))
    results_cross <- c(results_cross, as.vector(cross_dists))
    
    groups_inside <- c(groups_inside, rep(g, length(as.vector(inside_dists))))
    groups_cross <- c(groups_cross, rep(g, length(as.vector(cross_dists))))
    
  }
  
  inside_df <- data.frame(
    "group_name" = unlist(as.character(groups_inside)),
    "dist" = unlist(results_inside)
  )
  
  cross_df <- data.frame(
    "group_name" = unlist(as.character(groups_cross)),
    "dist" = unlist(results_cross)
  )
  
  inside_df$group_name <- as.factor(inside_df$group_name)
  cross_df$group_name <- as.factor(cross_df$group_name)
  
  results_df <- rbind(
    data.frame(
      "group_name"=inside_df$group_name, "dist"=inside_df$dist, "taxonomy"="same group"
    ),
    data.frame(
      "group_name"=cross_df$group_name, "dist"=cross_df$dist, "taxonomy"="different group"
    )
  )
  results_df$group_name <- as.factor(results_df$group_name)
  results_df$taxonomy <- as.factor(results_df$taxonomy)
  
  return(results_df)
  
}

boxplot_distances <- function(tax_results, avg_value, x_label, site, title="", include_not_in_dict=FALSE){
  
  if(site=="RA"){
    site_color <- "#F8766D"
  }else if(site=="RD"){
    site_color <- "#00BA38"
  }else if(site=="RL"){
    site_color <- "#619CFF"
  }
  
  if(x_label=="Phylum"){
    label_dict <- phylum_dict
  }else if(x_label=="Class"){
    label_dict <- class_dict
  }else if(x_label=="Order"){
    label_dict <- data.frame("ids"=c(), "names"=c())
  }else{
    label_dict <- data.frame("ids"=c(), "names"=c())
  }
  
  numerous_groups <- names(table(tax_results[tax_results$taxonomy=="same group", "group_name"]))[table(tax_results[tax_results$taxonomy=="same group","group_name"]) > 100]
  
  tax_results <- tax_results[tax_results$group_name %in% numerous_groups, ]
  tax_results <- tax_results[tax_results$group_name != -1, ]
  tax_results$group_name_ncbi <- ""
  for(g in unique(tax_results$group_name)){
    if(g %in% label_dict$ids){
      tax_results[tax_results$group_name==g, "group_name_ncbi"] <- label_dict[label_dict$ids==g, "names"]
    }else{
      message(g, " not in dictionary")
      if(include_not_in_dict){
        tax_results[tax_results$group_name==g, "group_name_ncbi"] <- as.character(g)
      }
    }
  }
  
  if(title!=""){
    title <- paste(site, " site - ", title, sep='')
  }else{
    title <- paste(site, " site", sep='')
  }
  
  p <- ggplot(
    tax_results[tax_results$group_name %in% numerous_groups, ]
  ) + geom_vline(
    aes(xintercept=avg_value), color="red", linewidth=1
  ) + geom_boxplot(
    aes(y=group_name_ncbi, x=dist, fill=taxonomy)
  ) + scale_discrete_manual(
    aesthetic="fill", values=c(site_color,"gray"), breaks=c("same group", "different group"), labels=c("same group", "different group")
  ) + labs(
    title=title,
    y=x_label,
    x="Distance"
  ) + theme_minimal(
    #  ) + theme(
    #    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )
  
  return(p)
  
}

#-------------------------------------------------------------------------------

RA_phylum_results <- get_distributions(RA_lcc_shortestPaths, "phylum_id")
RD_phylum_results <- get_distributions(RD_lcc_shortestPaths, "phylum_id")
RL_phylum_results <- get_distributions(RL_lcc_shortestPaths, "phylum_id")

RA_phylum_plot <- boxplot_distances(RA_phylum_results, RA_mean, "Phylum", "RA")
RD_phylum_plot <- boxplot_distances(RD_phylum_results, RD_mean, "Phylum", "RD")
RL_phylum_plot <- boxplot_distances(RL_phylum_results, RL_mean, "Phylum", "RL")

title <- "Shortest paths between species"

grid.arrange(
  top=title,
  RA_phylum_plot, RD_phylum_plot, RL_phylum_plot,
  ncol=1
)

#-------------------------------------------------------------------------------

RA_class_results <- get_distributions(RA_lcc_shortestPaths, "class_id")
RD_class_results <- get_distributions(RD_lcc_shortestPaths, "class_id")
RL_class_results <- get_distributions(RL_lcc_shortestPaths, "class_id")

RA_class_plot <- boxplot_distances(RA_class_results, RA_mean, "Class", "RA")
RD_class_plot <- boxplot_distances(RD_class_results, RD_mean, "Class", "RD")
RL_class_plot <- boxplot_distances(RL_class_results, RL_mean, "Class", "RL")

RA_order_results <- get_distributions(RA_lcc_shortestPaths, "order_id")
RD_order_results <- get_distributions(RD_lcc_shortestPaths, "order_id")
RL_order_results <- get_distributions(RL_lcc_shortestPaths, "order_id")

RA_order_plot <- boxplot_distances(RA_order_results, RA_mean, "Order", "RA", include_not_in_dict=TRUE)
RD_order_plot <- boxplot_distances(RD_order_results, RD_mean, "Order", "RD", include_not_in_dict=TRUE)
RL_order_plot <- boxplot_distances(RL_order_results, RL_mean, "Order", "RL", include_not_in_dict=TRUE)

title <- "Shortest paths between species"

grid.arrange(
  top=title,
  RA_class_plot, RA_order_plot,
  RD_class_plot, RD_order_plot,
  RL_class_plot, RL_order_plot,
  ncol=2
)

#-------------------------------------------------------------------------------

RA_phylum_results$site <- "RA"
RD_phylum_results$site <- "RD"
RL_phylum_results$site <- "RL"
RA_class_results$site <- "RA"
RD_class_results$site <- "RD"
RL_class_results$site <- "RL"
RA_order_results$site <- "RA"
RD_order_results$site <- "RD"
RL_order_results$site <- "RL"
RA_phylum_results$lvl <- "phylum"
RD_phylum_results$lvl <- "phylum"
RL_phylum_results$lvl <- "phylum"
RA_class_results$lvl <- "class"
RD_class_results$lvl <- "class"
RL_class_results$lvl <- "class"
RA_order_results$lvl <- "order"
RD_order_results$lvl <- "order"
RL_order_results$lvl <- "order"

all_results <- rbind(
  RA_phylum_results, RD_phylum_results, RL_phylum_results,
  RA_class_results, RD_class_results, RL_class_results,
  RA_order_results, RD_order_results, RL_order_results
)

all_results$lvl <- factor(all_results$lvl, levels=c("phylum", "class", "order"))

p_RA <- ggplot(
  all_results[all_results$site=="RA", ]
) + geom_vline(
  aes(xintercept=RA_mean), color="red", linewidth=1
) + geom_boxplot(
  aes(y=lvl, x=dist, fill=taxonomy)
) + scale_discrete_manual(
  aesthetic="fill", values=c("#F8766D","gray"), breaks=c("same group", "different group"), labels=c("same group", "different group")
) + labs(
  title="RA",
  y="Taxonomic level",
  x="Distance"
) + theme_minimal(
)

p_RD <- ggplot(
  all_results[all_results$site=="RD", ]
) + geom_vline(
  aes(xintercept=RD_mean), color="red", linewidth=1
) + geom_boxplot(
  aes(y=lvl, x=dist, fill=taxonomy)
) + scale_discrete_manual(
  aesthetic="fill", values=c("#00BA38","gray"), breaks=c("same group", "different group"), labels=c("same group", "different group")
) + labs(
  title="RD",
  y="Taxonomic level",
  x="Distance"
) + theme_minimal(
)

p_RL <- ggplot(
  all_results[all_results$site=="RL", ]
) + geom_vline(
  aes(xintercept=RL_mean), color="red", linewidth=1
) + geom_boxplot(
  aes(y=lvl, x=dist, fill=taxonomy)
) + scale_discrete_manual(
  aesthetic="fill", values=c("#619CFF","gray"), breaks=c("same group", "different group"), labels=c("same group", "different group")
) + labs(
  title="RL",
  y="Taxonomic level",
  x="Distance"
) + theme_minimal(
)

grid.arrange(
  top="Shortest paths between species, within same taxonomical groups and between different taxonomical groups",
  p_RA, p_RD, p_RL,
  ncol=1
)




