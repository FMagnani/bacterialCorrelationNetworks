library(data.table)
library(ggplot2)
library(VGAM)
library(igraph)

load(path_correlations)
globalThreshold_metrics <- as.data.frame(fread(path_globalThreshold_metrics))

RA_corr <- correlations$RA
RD_corr <- correlations$RD
RL_corr <- correlations$RL

#-------------------------------------------------------------------------------
# DENSITY & CORRELATION THRESHOLDS

as <- c(1.8, 1.9, 2.0, 2.1, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.6, 2.8, 3.0, 3.5)

zscore_to_density <- function(site_corr, as){
  ds <- c()
  for(a in as){
    with_filter <- c()
    for(r in rownames(site_corr)){
      this_row <- site_corr[r, ]
      m <- median(this_row)
      std <- sd(this_row)
      with_filter <- c(with_filter, sum(abs(this_row) > m+a*std))
    }
    ds <- c(ds, sum(100*sum(with_filter)/(1975*1975)))  
    message(a)
  }
  return(ds)
}

RA_densities <- zscore_to_density(RA_corr, as)
RD_densities <- zscore_to_density(RD_corr, as)
RL_densities <- zscore_to_density(RL_corr, as)

d_df <- rbind(
  data.frame("a"=unlist(as), "d"=unlist(RA_densities), "site"="RA"),
  data.frame("a"=unlist(as), "d"=unlist(RD_densities), "site"="RD"),
  data.frame("a"=unlist(as), "d"=unlist(RL_densities), "site"="RL")
)

p <- ggplot(
  d_df
) + geom_vline(
  aes(xintercept=1), color="red", lty="dashed", linewidth=1.5
) + geom_line(
  aes(x=d, y=a, group=site), colour="black", linewidth=1.5, 
) + geom_point(
  aes(x=d, y=a, fill=site), shape=21, colour="black", size=5
) + labs(
  title="Z-score threshold and resulting network density",
  x="Density",
  y="Z-score threshold"
) + theme_minimal()
p

#-------------------------------------------------------------------------------

make_adj <- function(site_corr, a){

  site_links <- matrix(0, 1975, 1975)
  rownames(site_links) <- rownames(site_corr)
  colnames(site_links) <- colnames(site_corr)
  for(r in rownames(site_links)){
    this_row <- site_corr[r, ]
    m <- median(this_row)
    std <- sd(this_row)
    pos_links <- (this_row > m+a*std)
    neg_links <- (this_row < m-a*std)
    site_links[r, pos_links] <- +1
    site_links[r, neg_links] <- -1
  }
  
  # Make symmetric
  mask_ones <- apply(site_links==t(site_links), c(1,2), as.integer)
  site_links <- site_links*mask_ones

  return(site_links)

}

# Save adjacencies to disk
save_dir <- path_localThresholdData
for(site_name in c("RA", "RD", "RL")){
  site_corr <- correlations[[site_name]]
  for(a in as){
    adj <- make_adj(site_corr, a)
    save_name <- paste(save_dir, site_name, "_", as.character(a), ".csv", sep='')
    write.csv(adj, save_name)
    message(save_name)
  }
}

read_adj <- function(site_name, a){
  file_name <- paste(save_dir, site_name, "_", as.character(a), ".csv", sep='')
  t <- as.data.frame(fread(file_name, header=TRUE))
  rownames(t) <- t$V1
  t$V1 <- NULL
  return(t)
}

#-------------------------------------------------------------------------------

RA_RD_cl <- c()
RA_RL_cl <- c()
RD_RL_cl <- c()
all_cl <- c()
for(a in as){
  
  RA_adj <- read_adj("RA", a)
  RD_adj <- read_adj("RD", a)
  RL_adj <- read_adj("RL", a)
  
  RA_Nlinks <- sum(abs(RA_adj))
  RD_Nlinks <- sum(abs(RD_adj))
  RL_Nlinks <- sum(abs(RL_adj))
  
  RA_RD_tot_links <- sum((abs(RA_adj)+abs(RD_adj))!=0)
  RA_RL_tot_links <- sum((abs(RA_adj)+abs(RL_adj))!=0)
  RD_RL_tot_links <- sum((abs(RD_adj)+abs(RL_adj))!=0)
  all_tot_links <- sum((abs(RA_adj)+abs(RD_adj)+abs(RL_adj))!=0)
  
  RA_adj[RA_adj==0] <- "RA"
  RD_adj[RD_adj==0] <- "RD"
  RL_adj[RL_adj==0] <- "RL"
  
  RA_RD_common_links <- 100*sum(RA_adj==RD_adj)/RA_RD_tot_links
  RA_RL_common_links <- 100*sum(RA_adj==RL_adj)/RA_RL_tot_links
  RD_RL_common_links <- 100*sum(RD_adj==RL_adj)/RD_RL_tot_links
  all_common_links <- 100*sum((RA_adj==RD_adj)&(RA_adj==RL_adj))/all_tot_links
  
  RA_RD_cl <- c(RA_RD_cl, RA_RD_common_links)
  RA_RL_cl <- c(RA_RL_cl, RA_RL_common_links)
  RD_RL_cl <- c(RD_RL_cl, RD_RL_common_links)
  all_cl <- c(all_cl, all_common_links)
  
  message(a)
  
}

localThreshold_metrics <- list("RA_RD_cl"=RA_RD_cl, "RA_RL_cl"=RA_RL_cl, "RD_RL_cl"=RD_RL_cl, "all_cl"=all_cl)
localThreshold_metrics$RA_density <- unlist(RA_densities)
localThreshold_metrics$RD_density <- unlist(RD_densities)
localThreshold_metrics$RL_density <- unlist(RL_densities)
localThreshold_metrics$RA_RD_density <- (localThreshold_metrics$RA_density+localThreshold_metrics$RD_density)/2
localThreshold_metrics$RA_RL_density <- (localThreshold_metrics$RA_density+localThreshold_metrics$RL_density)/2
localThreshold_metrics$RD_RL_density <- (localThreshold_metrics$RD_density+localThreshold_metrics$RL_density)/2
localThreshold_metrics$avg_density <- (localThreshold_metrics$RA_density+localThreshold_metrics$RD_density+localThreshold_metrics$RL_density)/3

#-------------------------------------------------------------------------------

cls_comparison_df <- rbind(
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cl"=unlist(globalThreshold_metrics$RA_RD_cl), "Thresholding"="global", "couple"="G_1"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cl"=unlist(globalThreshold_metrics$RA_RL_cl), "Thresholding"="global", "couple"="G_2"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cl"=unlist(globalThreshold_metrics$RD_RL_cl), "Thresholding"="global", "couple"="G_3"),
  data.frame("density"=unlist(localThreshold_metrics$RA_RD_density), "cl"=unlist(localThreshold_metrics$RA_RD_cl), "Thresholding"="local", "couple"="L_1"),
  data.frame("density"=unlist(localThreshold_metrics$RA_RL_density), "cl"=unlist(localThreshold_metrics$RA_RL_cl), "Thresholding"="local", "couple"="L_2"),
  data.frame("density"=unlist(localThreshold_metrics$RD_RL_density), "cl"=unlist(localThreshold_metrics$RD_RL_cl), "Thresholding"="local", "couple"="L_3")
)
cls_comparison_df$Thresholding <- as.factor(cls_comparison_df$Thresholding)
cls_comparison_df$couple <- as.factor(cls_comparison_df$couple)

p <- ggplot(
  cls_comparison_df[cls_comparison_df$density<3, ]
) + geom_line(
  aes(x=density, y=cl, colour=Thresholding, group=couple), linewidth=2, 
) + geom_point(
  aes(x=density, y=cl, fill=Thresholding), shape=21, colour="black", size=5
) + scale_discrete_manual(
  aesthetic="colour", values=c("navy","red4"), breaks=c("global", "local"), labels=c("global", "local")
) + scale_discrete_manual(
  aesthetic="fill", values=c("lightblue","red"), breaks=c("global", "local"), labels=c("global", "local")
) + labs(
  title="Percentage of common links after thresholding",
  x="Density",
  y="Percentage of common links"
) + theme_minimal()
p

#-------------------------------------------------------------------------------

lcc <- function(adj){
  
  g <- graph_from_adjacency_matrix(
    as.matrix(adj), 
    mode = "undirected", weighted = TRUE, diag = FALSE, 
    add.colnames = NA, add.rownames = NULL
  )
  
  comps <- components(g)
  lcc_index <- which(comps$csize==max(comps$csize))
  if(length(lcc_index)>1){
    message("Warning: More LCCs with same size are present!")
    lcc_index <- lcc_index[1]
  }
  nodes <- rownames(adj)
  lcc_nodes <- nodes[comps$membership==lcc_index]
  
  return(lcc_nodes)
  
}

RA_RD_corrDegree <- c()  
RA_RL_corrDegree <- c()  
RD_RL_corrDegree <- c() 
RA_slcc <- c()
RD_slcc <- c()
RL_slcc <- c()
RA_RD_clcc <- c()
RA_RL_clcc <- c()
RD_RL_clcc <- c()
all_clcc <- c()
for(a in as){
  
  RA_adj <- read_adj("RA", a)
  RD_adj <- read_adj("RD", a)
  RL_adj <- read_adj("RL", a)
  
  RA_degree <- rowSums(abs(RA_adj))
  RD_degree <- rowSums(abs(RD_adj))
  RL_degree <- rowSums(abs(RL_adj))
  
  # Correlation of degree
  RA_RD_corrDegree <- c(RA_RD_corrDegree, cor(RA_degree, RD_degree))  
  RA_RL_corrDegree <- c(RA_RL_corrDegree, cor(RA_degree, RL_degree))  
  RD_RL_corrDegree <- c(RD_RL_corrDegree, cor(RD_degree, RL_degree))  
  
  # LCC
  RA_lcc <- lcc(RA_adj)
  RD_lcc <- lcc(RD_adj)
  RL_lcc <- lcc(RL_adj)
  
  # Size LCC
  RA_size_lcc <- length(RA_lcc)
  RD_size_lcc <- length(RD_lcc)
  RL_size_lcc <- length(RL_lcc)
  
  RA_slcc <- c(RA_slcc, RA_size_lcc)
  RD_slcc <- c(RD_slcc, RD_size_lcc)
  RL_slcc <- c(RL_slcc, RL_size_lcc)
  
  # Common LCC
  RA_RD_common_lcc <- length(intersect(RA_lcc, RD_lcc))
  RA_RL_common_lcc <- length(intersect(RA_lcc, RL_lcc))
  RD_RL_common_lcc <- length(intersect(RD_lcc, RL_lcc))
  all_common_lcc <- length(intersect(intersect(RA_lcc, RD_lcc), RL_lcc))
  
  RA_RD_clcc <- c(RA_RD_clcc, RA_RD_common_lcc)
  RA_RL_clcc <- c(RA_RL_clcc, RA_RL_common_lcc)
  RD_RL_clcc <- c(RD_RL_clcc, RD_RL_common_lcc)
  all_clcc <- c(all_clcc, all_common_lcc)
  
  message(a)
  
}

localThreshold_metrics$RA_slcc <- unlist(RA_slcc)
localThreshold_metrics$RD_slcc <- unlist(RD_slcc)
localThreshold_metrics$RL_slcc <- unlist(RL_slcc)

localThreshold_metrics$RA_RD_corrDegree <- unlist(RA_RD_corrDegree)
localThreshold_metrics$RA_RL_corrDegree <- unlist(RA_RL_corrDegree)
localThreshold_metrics$RD_RL_corrDegree <- unlist(RD_RL_corrDegree)

localThreshold_metrics$RA_RD_clcc <- unlist(RA_RD_clcc)
localThreshold_metrics$RA_RL_clcc <- unlist(RA_RL_clcc)
localThreshold_metrics$RD_RL_clcc <- unlist(RD_RL_clcc)
localThreshold_metrics$all_clcc <- unlist(all_clcc)

#-------------------------------------------------------------------------------

# Correlation degree: Pearson
corrdegree_comparison_df <- rbind(
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cd"=unlist(globalThreshold_metrics$RA_RD_corrDegree), "Thresholding"="global", "couple"="G_RA-RD"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cd"=unlist(globalThreshold_metrics$RA_RL_corrDegree), "Thresholding"="global", "couple"="G_RA-RL"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cd"=unlist(globalThreshold_metrics$RD_RL_corrDegree), "Thresholding"="global", "couple"="G_RD-RL"),
  data.frame("density"=unlist(localThreshold_metrics$RA_RD_density), "cd"=unlist(localThreshold_metrics$RA_RD_corrDegree), "Thresholding"="local", "couple"="L_RA-RD"),
  data.frame("density"=unlist(localThreshold_metrics$RA_RL_density), "cd"=unlist(localThreshold_metrics$RA_RL_corrDegree), "Thresholding"="local", "couple"="L_RA-RL"),
  data.frame("density"=unlist(localThreshold_metrics$RD_RL_density), "cd"=unlist(localThreshold_metrics$RD_RL_corrDegree), "Thresholding"="local", "couple"="L_RD-RL")
)

p <- ggplot(
  corrdegree_comparison_df[corrdegree_comparison_df$density>0.05, ]
) + geom_line(
  aes(x=density, y=cd, colour=Thresholding, group=couple), linewidth=2, 
) + geom_point(
  aes(x=density, y=cd, fill=Thresholding), shape=21, colour="black", size=5
) + scale_discrete_manual(
  aesthetic="colour", values=c("navy","red4"), breaks=c("global", "local"), labels=c("global", "local")
) + scale_discrete_manual(
  aesthetic="fill", values=c("lightblue","red"), breaks=c("global", "local"), labels=c("global", "local")
) + labs(
  title="Pearson correlation of degrees of nodes after thresholding",
  x="Density",
  y="Correlations"
) + theme_minimal()
p

# Size LCC
sizeLCC_comparison_df <- rbind(
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "sLCC"=unlist(globalThreshold_metrics$RA_slcc), "thr"="global", "site"="RA", "group"=1),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "sLCC"=unlist(globalThreshold_metrics$RD_slcc), "thr"="global", "site"="RD", "group"=2),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "sLCC"=unlist(globalThreshold_metrics$RL_slcc), "thr"="global", "site"="RL", "group"=3),
  data.frame("density"=unlist(localThreshold_metrics$RA_density), "sLCC"=unlist(localThreshold_metrics$RA_slcc), "thr"="local", "site"="RA", "group"=4),
  data.frame("density"=unlist(localThreshold_metrics$RD_density), "sLCC"=unlist(localThreshold_metrics$RD_slcc), "thr"="local", "site"="RD", "group"=5),
  data.frame("density"=unlist(localThreshold_metrics$RL_density), "sLCC"=unlist(localThreshold_metrics$RL_slcc), "thr"="local", "site"="RL", "group"=6)
)

p <- ggplot(
  sizeLCC_comparison_df[sizeLCC_comparison_df$density<3, ]
) + geom_line(
  aes(x=density, y=sLCC, colour=thr, group=group), linewidth=2, 
) + geom_point(
  aes(x=density, y=sLCC, fill=site), shape=21, colour="black", size=5
) + scale_discrete_manual(
  aesthetic="colour", values=c("navy","red4"), breaks=c("global", "local"), labels=c("global", "local")
) + scale_discrete_manual(
  aesthetic="fill", values=c("#F8766D","#00BA38", "#619CFF"), breaks=c("RA", "RD", "RL"), labels=c("RA", "RD", "RL")
) + labs(
  title="Size of Largest Connected Component",
  x="Density",
  y="Number of nodes"
) + theme_minimal()
p

# Common LCC
corrdegree_comparison_df <- rbind(
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "clcc"=unlist(globalThreshold_metrics$RA_RD_clcc), "Thresholding"="global", "couple"="G_RA-RD"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "clcc"=unlist(globalThreshold_metrics$RA_RL_clcc), "Thresholding"="global", "couple"="G_RA-RL"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "clcc"=unlist(globalThreshold_metrics$RD_RL_clcc), "Thresholding"="global", "couple"="G_RD-RL"),
  data.frame("density"=unlist(localThreshold_metrics$RA_RD_density), "clcc"=unlist(localThreshold_metrics$RA_RD_clcc), "Thresholding"="local", "couple"="L_RA-RD"),
  data.frame("density"=unlist(localThreshold_metrics$RA_RL_density), "clcc"=unlist(localThreshold_metrics$RA_RL_clcc), "Thresholding"="local", "couple"="L_RA-RL"),
  data.frame("density"=unlist(localThreshold_metrics$RD_RL_density), "clcc"=unlist(localThreshold_metrics$RD_RL_clcc), "Thresholding"="local", "couple"="L_RD-RL")
)

p <- ggplot(
  corrdegree_comparison_df
) + geom_line(
  aes(x=density, y=clcc, colour=Thresholding, group=couple), linewidth=2, 
) + geom_point(
  aes(x=density, y=clcc, fill=Thresholding), shape=21, colour="black", size=5
) + scale_discrete_manual(
  aesthetic="colour", values=c("navy","red4"), breaks=c("global", "local"), labels=c("global", "local")
) + scale_discrete_manual(
  aesthetic="fill", values=c("lightblue","red"), breaks=c("global", "local"), labels=c("global", "local")
) + labs(
  title="Pearson correlation of degrees of nodes after thresholding",
  x="Density",
  y="Correlations"
) + theme_minimal()
p

#-------------------------------------------------------------------------------

write.csv(localThreshold_metrics, path_localThreshold_metrics)
