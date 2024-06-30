library(igraph)
library(mgnet) # https://github.com/Fuschi/mgnet
library(stringr)

load(path_adjacencies)

get_communities <- function(adj, resistance){
  
  g <- graph_from_adjacency_matrix(adj, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA, add.rownames = NULL)
  communities <- custom_cluster_signed(g)
  
  return(communities)
  
}

custom_cluster_signed <- function(obj, Resistance=0, Penalty_Coefficien=1, add.names=TRUE, routine="l 1"){  
  
  OS <- "Linux"
  graph <- obj
  
  #Check Graph
  if(!is.weighted(graph)) stop("graph must be weighted graph")
  if(is.directed(graph))  stop("graph must be undirected")
  if(file.exists("graph.net")) file.remove("graph.net")
  if(!is.numeric(Resistance)) stop("Resistance must be numeric")
  if(!is.numeric(Penalty_Coefficien)) stop("Penalty_Coefficien must be a number >= 0")
  if(Penalty_Coefficien<0) stop("Penalty_Coefficien must be a number >= 0")
  if(add.names & !igraph::is_named(obj)) stop("graph has not vertices names attribute")
  
  adj <- as_adjacency_matrix(graph, attr="weight", sparse=FALSE)

  #get path of executable Communities_Detection.exe
  path <- system.file("exec", package="mgnet", mustWork=TRUE)
  path <- paste(path,"/",sep="")
  
  #write graph in pajek format as request from executable
  write_graph(graph  = graph_from_adjacency_matrix(adjmatrix=adj, mode="undirected", weighted=TRUE),
              file   = paste(path,"graph.net",sep=""),
              format ="pajek")
  
  #path for graph and results files
  path.graph <- paste(path,"graph.net",sep="")
  path.result <- paste(path, "res.txt", sep="")
  
  #save command for the execution
  cmd <- paste(path,"Communities_Detection_Linux.exe N WS ", routine, sep="")
  cmd <- paste(cmd, Resistance, Penalty_Coefficien, path.graph, path.result, sep=" ")
  
  system(cmd,ignore.stdout=TRUE,ignore.stderr=TRUE)
  
  #open file of results
  res.file = file(path.result, "r")
  
  #read and store results
  file.lines <- readLines(res.file)
  info <- file.lines[2]
  modularity <- as.numeric(str_split(file.lines[3]," ")[[1]][3])
  
  vertex.num <- as.numeric(str_split(file.lines[5]," ")[[1]][4])
  comm.num <- as.numeric(str_split(file.lines[6]," ")[[1]][4])
  comm.vert <- vector(mode="list",length=comm.num)
  
  for(line in 8:(7+comm.num)){
    
    i <- line-7
    #comm[[i]]
    tmp <- unlist(str_split(file.lines[line]," "))
    tmp <- as.numeric(tmp[2:length(tmp)])
    
    comm.vert[[i]] <- tmp
  }
  
  comm <- vector(mode="integer", length=vertex.num)
  for(c in 1:comm.num){comm[comm.vert[[c]]] <- c}
  
  #close connection
  close(res.file)
  
  #remove files
  file.remove(path.graph)
  file.remove(path.result)
  
  res <- make_clusters(graph=graph,
                       membership=comm,
                       algorithm="signed weights louvain",
                       modularity=modularity)
  
  isolated.membership <- which(sizes(res)==1)
  res$membership[res$membership %in% isolated.membership] <- 0
  
  if(add.names) names(res$membership) <- V(graph)$name
  
  #return the results as communities structure of igraph package
  return(res)
}

#---------------------------------------------------------------------------------------------

resistance <- 0

RA_g <- graph_from_adjacency_matrix(adjacencies$RA_adj, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA, add.rownames = NULL)
RD_g <- graph_from_adjacency_matrix(adjacencies$RD_adj, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA, add.rownames = NULL)
RL_g <- graph_from_adjacency_matrix(adjacencies$RL_adj, mode = "undirected", weighted = TRUE, diag = FALSE, add.colnames = NA, add.rownames = NULL)

# routine="r-s-e-ll!rfr-trfr 5"
RA_comm <- custom_cluster_signed(RA_g)
RD_comm <- custom_cluster_signed(RD_g)
RL_comm <- custom_cluster_signed(RL_g)

#---------------------------------------------------------------------------------------------

# For UMAP (py)
write.csv(RA_comm$membership, "RA_communities.csv")
write.csv(RD_comm$membership, "RD_communities.csv")
write.csv(RL_comm$membership, "RL_communities.csv")

#---------------------------------------------------------------------------------------------

# For R 
communities <- list(
  "RL"=RL_comm$membership, "RA"=RA_comm$membership, "RD"=RD_comm$membership
)
save(communities, file="communities.RData")
