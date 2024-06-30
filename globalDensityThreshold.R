library(data.table)
library(ggplot2)
library(VGAM)
library(igraph)

load(path_correlations)

RA_corr <- correlations$RA
RD_corr <- correlations$RD
RL_corr <- correlations$RL

#-------------------------------------------------------------------------------
# DENSITY & CORRELATION THRESHOLDS

densities <- seq(2, 50, 2)/1000

density_to_correlation <- function(site_corr, densities){
  corr_thresholds <- c()
  for(d in densities){
    N_links <- length(site_corr)
    index_thr <- 1 + as.integer(d*N_links)
    corr_thr <- abs(
      site_corr[
        order(abs(
          site_corr
        ), decreasing=TRUE)
      ][[index_thr]]
    )
    corr_thresholds <- c(corr_thresholds, corr_thr)
  }
  return(corr_thresholds)
}

RA_corr_thrs <- density_to_correlation(RA_corr, densities)
RD_corr_thrs <- density_to_correlation(RD_corr, densities)
RL_corr_thrs <- density_to_correlation(RL_corr, densities)

d_df <- rbind(
  data.frame("d"=100*unlist(densities), "c"=unlist(RA_corr_thrs), "site"="RA"),
  data.frame("d"=100*unlist(densities), "c"=unlist(RD_corr_thrs), "site"="RD"),
  data.frame("d"=100*unlist(densities), "c"=unlist(RL_corr_thrs), "site"="RL")
)

p <- ggplot(
  d_df[d_df$c!=0, ]
) + geom_vline(
  aes(xintercept=1), color="red", lty="dashed", linewidth=1.5
) + geom_line(
  aes(x=d, y=c, group=site), colour="black", linewidth=2, 
) + geom_point(
  aes(x=d, y=c, fill=site), shape=21, colour="black", size=5
) + labs(
  title="Correlation threshold and resulting network density",
  x="Density",
  y="Correlation threshold"
) + theme_minimal()
p

#-------------------------------------------------------------------------------

make_adj <- function(site_corr, corr_thr){
  adj <- site_corr
  adj[abs(adj)<corr_thr] <- 0
  adj[adj > 0] <-  1
  adj[adj < 0] <- -1
  return(adj)
}

RA_RD_cl <- c()
RA_RL_cl <- c()
RD_RL_cl <- c()
all_cl <- c()
for(idx in 1:length(RA_corr_thrs)){
  
  RA_adj <- make_adj(RA_corr, RA_corr_thrs[idx])
  RD_adj <- make_adj(RD_corr, RD_corr_thrs[idx])
  RL_adj <- make_adj(RL_corr, RL_corr_thrs[idx])
  
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
  
  message(idx, "    ", densities[[idx]])

}

globalThreshold_metrics <- list("RA_RD_cl"=RA_RD_cl, "RA_RL_cl"=RA_RL_cl, "RD_RL_cl"=RD_RL_cl, "all_cl"=all_cl)
globalThreshold_metrics$density <- unlist(densities)

#-------------------------------------------------------------------------------

cls_df <- rbind(
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cl"=unlist(globalThreshold_metrics$RA_RD_cl), "group"="RA-RD", "comparison"="couple"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cl"=unlist(globalThreshold_metrics$RA_RL_cl), "group"="RA-RL", "comparison"="couple"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cl"=unlist(globalThreshold_metrics$RD_RL_cl), "group"="RD-RL", "comparison"="couple"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cl"=unlist(globalThreshold_metrics$all_cl), "group"="all", "comparison"="all")
)

p <- ggplot(
  cls_df
) + annotate(
  "rect", 
  colour="red", fill="grey", alpha=0.6, 
  xmin = 0.35, xmax = 1.15,
  #ymin = 20, ymax = 23.5
  ymin = -Inf, ymax = Inf
) + geom_line(
  aes(x=density, y=cl, group=group, colour=group), linewidth=2, 
) + geom_point(
  aes(x=density, y=cl, fill=group), shape=21, size=5
) + scale_colour_manual(
  values = c("black", "#805B56", "#8077D5", "#00D180")
) + scale_fill_manual(
  values = c("black", "#805B56", "#8077D5", "#00D180")
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
for(idx in 1:length(RA_corr_thrs)){
  
  RA_adj <- make_adj(RA_corr, RA_corr_thrs[idx])
  RD_adj <- make_adj(RD_corr, RD_corr_thrs[idx])
  RL_adj <- make_adj(RL_corr, RL_corr_thrs[idx])
  
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
  
  message(idx, "    ", densities[[idx]])
  
}

globalThreshold_metrics$RA_slcc <- unlist(RA_slcc)
globalThreshold_metrics$RD_slcc <- unlist(RD_slcc)
globalThreshold_metrics$RL_slcc <- unlist(RL_slcc)

globalThreshold_metrics$RA_RD_corrDegree <- unlist(RA_RD_corrDegree)
globalThreshold_metrics$RA_RL_corrDegree <- unlist(RA_RL_corrDegree)
globalThreshold_metrics$RD_RL_corrDegree <- unlist(RD_RL_corrDegree)

globalThreshold_metrics$RA_RD_clcc <- unlist(RA_RD_clcc)
globalThreshold_metrics$RA_RL_clcc <- unlist(RA_RL_clcc)
globalThreshold_metrics$RD_RL_clcc <- unlist(RD_RL_clcc)
globalThreshold_metrics$all_clcc <- unlist(all_clcc)

#-------------------------------------------------------------------------------

# Correlation degree: Pearson
corrdegree_df <- rbind(
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cd"=unlist(globalThreshold_metrics$RA_RD_corrDegree), "comparison"="RA-RD"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cd"=unlist(globalThreshold_metrics$RA_RL_corrDegree), "comparison"="RA-RL"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "cd"=unlist(globalThreshold_metrics$RD_RL_corrDegree), "comparison"="RD-RL")
)

p <- ggplot(
  corrdegree_df
) + geom_vline(
  aes(xintercept=1), color="red", lty="dashed", linewidth=1.5
) + geom_line(
  aes(x=density, y=cd, colour=comparison), linewidth=2, 
) + geom_point(
  aes(x=density, y=cd, fill=comparison), shape=21, colour="black", size=5
) + scale_colour_manual(
  values = c("#805B56", "#8077D5", "#00D180")
) + scale_fill_manual(
  values = c("#805B56", "#8077D5", "#00D180")
) + labs(
  title="Pearson correlation of nodes degree between couples of site-specific networks",
  x="Density",
  y="Correlation of nodes degree"
) + theme_minimal()
p

# Size LCC
sizeLCC_df <- rbind(
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "slcc"=unlist(globalThreshold_metrics$RA_slcc), "site"="RA"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "slcc"=unlist(globalThreshold_metrics$RD_slcc), "site"="RD"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "slcc"=unlist(globalThreshold_metrics$RL_slcc), "site"="RL")
)

p <- ggplot(
  sizeLCC_df
) + geom_vline(
  aes(xintercept=1), color="red", lty="dashed", linewidth=1.5
) + geom_line(
  aes(x=density, y=slcc, colour=site), linewidth=1.5, 
) + geom_point(
  aes(x=density, y=slcc, fill=site), shape=21, colour="black", size=3
) + labs(
  title="Size of the largest connected component of the site-specific networks",
  x="Density",
  y="Size of LCC (number of nodes)"
) + theme_minimal()
p

# Common LCC
commonLCC_df <- rbind(
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "clcc"=unlist(globalThreshold_metrics$RA_RD_clcc), "comparison"="RA-RD"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "clcc"=unlist(globalThreshold_metrics$RA_RL_clcc), "comparison"="RA-RL"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "clcc"=unlist(globalThreshold_metrics$RD_RL_clcc), "comparison"="RD-RL"),
  data.frame("density"=100*unlist(globalThreshold_metrics$density), "clcc"=unlist(globalThreshold_metrics$all_clcc), "comparison"="all")
)

p <- ggplot(
  commonLCC_df
) + geom_vline(
  aes(xintercept=1), color="red", lty="dashed", linewidth=1.5
) + geom_line(
  aes(x=density, y=clcc, colour=comparison), linewidth=2, 
) + geom_point(
  aes(x=density, y=clcc, fill=comparison), shape=21, colour="black", size=5
) + scale_colour_manual(
  values = c("black", "#805B56", "#8077D5", "#00D180")
) + scale_fill_manual(
  values = c("black", "#805B56", "#8077D5", "#00D180")
) + labs(
  title="Nodes of the LCC common to the site-specific networks",
  x="Density",
  y="Common nodes in LCC (number of nodes)"
) + theme_minimal()
p

#-------------------------------------------------------------------------------

write.csv(globalThreshold_metrics, path_globalThreshold_metrics)
