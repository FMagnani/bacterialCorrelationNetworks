library(igraph) 
library(mgnet) # https://github.com/Fuschi/mgnet
library(stringr)
library(data.table)

load(path_adjacencies)

taxtab <- read.csv(path_taxonomy_table)
taxtab$X <- NULL

load(path_communities)

#---------------------------------------------------------------------------------------------------------

species_csv <- as.data.frame(fread(path_species))
rownames_ids <- c()
for(n in rownames(adjacencies$RA_adj)){
  id <- species_csv[species_csv$names==n, "ncbi_ids"]
  rownames_ids <- c(rownames_ids, id)
}
rownames(adjacencies$RA_adj) <- rownames_ids
rownames(adjacencies$RD_adj) <- rownames_ids
rownames(adjacencies$RL_adj) <- rownames_ids
colnames(adjacencies$RA_adj) <- rownames_ids
colnames(adjacencies$RD_adj) <- rownames_ids
colnames(adjacencies$RL_adj) <- rownames_ids

#---------------------------------------------------------------------------------------------------------

taxonomic_modularity <- function(g, partition){
  
  OS <- "Linux"
  
  adj <- as_adjacency_matrix(g, attr="weight", sparse=FALSE)
  
  #get path of executable Communities_Detection.exe
  path <- system.file("exec", package="mgnet", mustWork=TRUE)
  path <- paste(path,"/",sep="")
  
  #path for graph and results files
  path.graph <- paste(path,"graph.net",sep="")
  path.partition <- paste(path,"partition.clu",sep="")
  path.result <- paste(path, "res.txt", sep="")
  
  #write graph in pajek format as request from executable
  write_graph(
    graph  = graph_from_adjacency_matrix(adjmatrix=adj, mode="undirected", weighted=TRUE),
    file   = path.graph,
    format ="pajek"
  )
  
  #write partition in pajek format as request from executable
  partition <- c(
    paste("*Vertices ", as.character(length(partition))), 
    partition
  )
  write.table(partition, path.partition, sep="", col.names=FALSE, row.names=FALSE, quote=FALSE)
  
  #save command for the execution
  cmd <- paste(path,"Modularity_Calculation.exe ", path.graph, " ", path.partition, " WS TC > ", path.result, sep='')
  #cmd <- paste(cmd, Resistance, Penalty_Coefficien, path.graph, path.result, sep=" ")
  
  system(cmd, ignore.stdout=FALSE, ignore.stderr=FALSE)
  
  #read and store results
  res <- read.table(path.result, header=F, sep="\t", fileEncoding="latin1")
  
  total_mod_line <- res[11, ]
  classes_mod_lines <- res[12:dim(res)[[1]], ]
  modularity_from_line <- function(line){return(as.numeric(strsplit(line, '=')[[1]][[2]]))}
  number_nodes_from_line <- function(line){
    return(
      as.numeric(strsplit(strsplit(strsplit(line, '=')[[1]][[3]], '[(]')[[1]][[2]], " ")[[1]][[1]])
    )
  }
  
  total_mod <- modularity_from_line(total_mod_line)
  classes_mod <- unlist(lapply(classes_mod_lines, modularity_from_line))
  classes_nodes <- unlist(lapply(classes_mod_lines, number_nodes_from_line))
  
  res <- data.frame(
    "tax"=unique(partition[2:length(partition)]),
    "modularity"=classes_mod,
    "numerosity"=classes_nodes
  )
  
  stopifnot(abs(total_mod - sum(res$modularity)) < 1e-04)
  node_numerosities <- table(partition)
  for(i in length(res$tax)){
    tax_id <- res[i, "tax"]
    num <- res[i, "numerosity"]
    stopifnot(num == node_numerosities[tax_id])
  }
  
  #remove files
  file.remove(path.graph)
  file.remove(path.partition)
  file.remove(path.result)
  
  return(res)
  
}

#---------------------------------------------------------------------------------------------------------

taxonomy_modularity <- list()
for(site in c("RA_adj", "RD_adj", "RL_adj")){
  site_results <- list()
  for(taxlvl in c("phylum", "class", "order", "family", "genus")){
    taxonomic_level <- paste(taxlvl, "_id", sep='')
    
    adj <- adjacencies[[site]]
    common_species <- intersect(rownames(adj), taxtab$species)
    
    common_taxtab <- taxtab[taxtab$species %in% common_species, ]
    adj <- adj[common_species, common_species]
    
    g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA, add.rownames = NULL)
    partition <- common_taxtab[[taxonomic_level]]
    
    res <- taxonomic_modularity(g, partition)
    
    site_results[[taxlvl]] <- res
    
    message(site, "    ", taxlvl)
  }
  taxonomy_modularity[[site]] <- site_results
}

#------------------------------------------------------------------------------------------------------------------

# Load adjacencies again since now we need names and not ncbi ids
load(path_adjacencies)

for(site in c("RA", "RD", "RL")){
  
  adj <- adjacencies[[site]]
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA, add.rownames = NULL)
  partition <- communities[[site]]
  
  res <- taxonomic_modularity(g, partition)
  
  taxonomy_modularity[[paste(site, "adj", sep='_')]]$net_community <- res
  
}

save(taxonomy_modularity, file=path_taxonomy_modularity)
