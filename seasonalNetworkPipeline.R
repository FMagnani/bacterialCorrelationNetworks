library(data.table)
library(matrixStats)
library(readr)
library(psych)
library(ggplot2)
library(gridExtra)
library(ggpubr)
library(VennDiagram)

load(path_bacteria_species_counts)

#-------------------------------------------------------------------------------

# FILTRAGGIO

get_row_stats <- function(df, site){
  nonzero_medians <- c()
  prevalence_and_numerosity <- c()
  for(n in rownames(df)){
    
    row <-  df[n, ]
    
    prevalence_and_numerosity <- c(
      prevalence_and_numerosity, 
      length(row[row >= 10])/length(row)
    )
    nonzero_filter <- row!=0
    if(length(row[nonzero_filter])!=0){
      nonzero_medians <- c(nonzero_medians, median(row[nonzero_filter]))
    }else{
      nonzero_medians <- c(nonzero_medians, 0)
    }
    
  }
  
  return(
    data.frame(
      "name" = unlist(rownames(df)),
      "nonzero_medians" = nonzero_medians,
      "prevalence_and_numerosity" = prevalence_and_numerosity,
      "site" = site
    )
  )
  
}

row_stats <- rbind(
  get_row_stats(bacteria_species_counts$RL, "RL"),
  get_row_stats(bacteria_species_counts$RA, "RA"),
  get_row_stats(bacteria_species_counts$RD, "RD")
)

#-------------------------------------------------------------------------------

# FILTERS

# Threshold on prevalence/numerosity
RL_prevnum_filter <- row_stats[row_stats$site=="RL", ]$prevalence_and_numerosity > 0.25
RA_prevnum_filter <- row_stats[row_stats$site=="RA", ]$prevalence_and_numerosity > 0.25
RD_prevnum_filter <- row_stats[row_stats$site=="RD", ]$prevalence_and_numerosity > 0.25

# Threshold on median
RL_median_filter <- row_stats[row_stats$site=="RL", ]$nonzero_medians > 5
RA_median_filter <- row_stats[row_stats$site=="RA", ]$nonzero_medians > 5
RD_median_filter <- row_stats[row_stats$site=="RD", ]$nonzero_medians > 5

ids <- rownames(bacteria_species_counts$RL)
RL_ids <- ids[RL_prevnum_filter & RL_median_filter]
RA_ids <- ids[RA_prevnum_filter & RA_median_filter]
RD_ids <- ids[RD_prevnum_filter & RD_median_filter]

common_bacteria <- intersect(intersect(RL_ids, RA_ids), RD_ids)

#-------------------------------------------------------------------------------

# CORREZIONE ZERI E CLR

correct_dataset <- function(site_counts){
  
  correct_counts <- site_counts
  for(date in names(site_counts)){
    
    sample <- site_counts[date]
    
    K  = length(sample)
    Z = length(sample[sample==0])
    
    delta <- 1/(K^2)
    
    nonzero_correction <- 1 - Z*delta/sum(sample)
    
    sample[sample==0] <- delta
    sample[sample!=0] <- nonzero_correction*sample[sample!=0]
    
    correct_counts[date] <- sample        
    
  }
  
  return(correct_counts)  
  
}

RL_filter <- bacteria_species_counts$RL[RL_ids, ]
RA_filter <- bacteria_species_counts$RA[RA_ids, ]
RD_filter <- bacteria_species_counts$RD[RD_ids, ]

RL_correct <- correct_dataset(RL_filter)
RA_correct <- correct_dataset(RA_filter)
RD_correct <- correct_dataset(RD_filter)

apply_CLR_to_dataset <- function(counts){
  # The dates should be in the columns
  # Each row should be a species/class/etc
  
  ref <- apply(counts, 2, function(x) mean(log(x)))
  transformed <- sweep(log(counts), 2, ref)
  #transformed <- as.matrix(log(counts) - ref)
  return(transformed)
  
}

RL_clr <- apply_CLR_to_dataset(RL_correct)
RA_clr <- apply_CLR_to_dataset(RA_correct)
RD_clr <- apply_CLR_to_dataset(RD_correct)

#-------------------------------------------------------------------------------

split_by_period <- function(site_clr, start_date, end_date){
  
  site_dates <- as.Date(names(site_clr), "%Y.%m.%d")
  
  site_clr <- site_clr[
    , site_dates>=start_date & site_dates<end_date
  ]
  
  return(site_clr)
  
}

cold_start <- as.Date("28.10.2019", "%d.%m.%Y")
cold_end <- as.Date("25.04.2020", "%d.%m.%Y")
hot_start <- as.Date("25.04.2020", "%d.%m.%Y")
hot_end <- as.Date("01.10.2020", "%d.%m.%Y")

RL_cold <- split_by_period(RL_clr, cold_start, cold_end)
RA_cold <- split_by_period(RA_clr, cold_start, cold_end)
RD_cold <- split_by_period(RD_clr, cold_start, cold_end)

RL_hot <- split_by_period(RL_clr, hot_start, hot_end)
RA_hot <- split_by_period(RA_clr, hot_start, hot_end)
RD_hot <- split_by_period(RD_clr, hot_start, hot_end)

#-------------------------------------------------------------------------------

get_corr_test <- function(site_data){
  test <- site_corr_test <- corr.test(
    x=t(site_data), 
    use="pairwise", method="pearson",
    adjust="bonferroni",
    ci=FALSE
  )
  return(test)
}

RA_cold_test <- get_corr_test(RA_cold)
RD_cold_test <- get_corr_test(RD_cold)
RL_cold_test <- get_corr_test(RL_cold)

RA_hot_test <- get_corr_test(RA_hot)
RD_hot_test <- get_corr_test(RD_hot)
RL_hot_test <- get_corr_test(RL_hot)

get_corr <- function(corr_test){
  site_corr <- corr_test$r
  diag(site_corr) <- 0
  return(site_corr)
}

RA_cold_corr <- get_corr(RA_cold_test)
RD_cold_corr <- get_corr(RD_cold_test)
RL_cold_corr <- get_corr(RL_cold_test)

RA_hot_corr <- get_corr(RA_hot_test)
RD_hot_corr <- get_corr(RD_hot_test)
RL_hot_corr <- get_corr(RL_hot_test)

RA_cold_corr <- RA_cold_corr[common_bacteria, common_bacteria]
RD_cold_corr <- RD_cold_corr[common_bacteria, common_bacteria]
RL_cold_corr <- RL_cold_corr[common_bacteria, common_bacteria]
RA_hot_corr <- RA_hot_corr[common_bacteria, common_bacteria]
RD_hot_corr <- RD_hot_corr[common_bacteria, common_bacteria]
RL_hot_corr <- RL_hot_corr[common_bacteria, common_bacteria]

#-------------------------------------------------------------------------------

cold_corr_df <- rbind(
  data.frame(
    "corr" = as.vector(unname(as.matrix(RA_cold_corr[upper.tri(RA_cold_corr)]))),
    "site" = "RA"
  ),
  data.frame(
    "corr" = as.vector(unname(as.matrix(RD_cold_corr[upper.tri(RD_cold_corr)]))),
    "site" = "RD"
  ),
  data.frame(
    "corr" = as.vector(unname(as.matrix(RL_cold_corr[upper.tri(RL_cold_corr)]))),
    "site" = "RL"
  )
)

hot_corr_df <- rbind(
  data.frame(
    "corr" = as.vector(unname(as.matrix(RA_hot_corr[upper.tri(RA_hot_corr)]))),
    "site" = "RA"
  ),
  data.frame(
    "corr" = as.vector(unname(as.matrix(RD_hot_corr[upper.tri(RD_hot_corr)]))),
    "site" = "RD"
  ),
  data.frame(
    "corr" = as.vector(unname(as.matrix(RL_hot_corr[upper.tri(RL_hot_corr)]))),
    "site" = "RL"
  )
)

# Correlation of correlations
round(
  cor(
    cold_corr_df[cold_corr_df$site=="RA", "corr"], 
    cold_corr_df[cold_corr_df$site=="RD", "corr"]
  ), digits=2
)

round(
  cor(
    cold_corr_df[cold_corr_df$site=="RA", "corr"], 
    cold_corr_df[cold_corr_df$site=="RL", "corr"]
  ), digits=2
)

round(
  cor(
    cold_corr_df[cold_corr_df$site=="RD", "corr"], 
    cold_corr_df[cold_corr_df$site=="RL", "corr"]
  ), digits=2
)

# HOT HOT
# RA-RD 0.55
# RA-RL 0.49 
# RD-RL 0.57

# HOT COLD
# RA-RD 0.33
# RA-RL 0.32
# RD-RL 0.35
# RA-RA 0.28
# RD-RD 0.35
# RL-RL 0.43

# COLD COLD
# RA-RD 0.64
# RA-RL 0.59
# RD-RL 0.64

#-------------------------------------------------------------------------------

density_to_correlation_threshold <- function(site_corr, density_thr){
  
  N_links <- length(site_corr)
  index_thr <- 1 + as.integer(density_thr*N_links)
  
  corr_thr <- abs(
    site_corr[
      order(abs(
        site_corr
      ), decreasing=TRUE)
    ][[index_thr]]
  )
  
  return(corr_thr)
  
}

RA_cold_r_thr <- density_to_correlation_threshold(RA_cold_corr, 1/100) # 0.90
RD_cold_r_thr <- density_to_correlation_threshold(RD_cold_corr, 1/100) # 0.89
RL_cold_r_thr <- density_to_correlation_threshold(RL_cold_corr, 1/100) # 0.85

RA_hot_r_thr <- density_to_correlation_threshold(RA_hot_corr, 1/100) # 0.84
RD_hot_r_thr <- density_to_correlation_threshold(RD_hot_corr, 1/100) # 0.84
RL_hot_r_thr <- density_to_correlation_threshold(RL_hot_corr, 1/100) # 0.82

make_adj <- function(site_corr, site_thr){
  site_adj <- site_corr
  site_adj[(abs(site_corr) < site_thr)] <- 0
  site_adj[site_adj > 0] <-  1
  site_adj[site_adj < 0] <- -1
  return(site_adj)
}

RA_cold_adj <- make_adj(RA_cold_corr, RA_cold_r_thr)
RD_cold_adj <- make_adj(RD_cold_corr, RD_cold_r_thr)
RL_cold_adj <- make_adj(RL_cold_corr, RL_cold_r_thr)

RA_hot_adj <- make_adj(RA_hot_corr, RA_hot_r_thr)
RD_hot_adj <- make_adj(RD_hot_corr, RD_hot_r_thr)
RL_hot_adj <- make_adj(RL_hot_corr, RL_hot_r_thr)

#-------------------------------------------------------------------------------

seasonal_adjacencies <- list(
  "cold"=list("RA"=RA_cold_adj, "RD"=RD_cold_adj, "RL"=RL_cold_adj),
  "hot"=list("RA"=RA_hot_adj, "RD"=RD_hot_adj, "RL"=RL_hot_adj)
)
save(seasonal_adjacencies, file="C:\\Users\\fede\\Desktop\\VEO\\DATA\\VEO_data\\seasonal_adjacencies.RData")

#-------------------------------------------------------------------------------

N_links <- sum(abs(RA_cold_adj))

both_one <- sum((seasonal_adjacencies$cold$RA==1)&(seasonal_adjacencies$cold$RD==1))
both_minus_one <- sum((seasonal_adjacencies$cold$RA==-1)&(seasonal_adjacencies$cold$RD==-1))
round(100*(both_one+both_minus_one)/N_links, digits=0)

both_one <- sum((seasonal_adjacencies$cold$RA==1)&(seasonal_adjacencies$cold$RL==1))
both_minus_one <- sum((seasonal_adjacencies$cold$RA==-1)&(seasonal_adjacencies$cold$RL==-1))
round(100*(both_one+both_minus_one)/N_links, digits=0)

both_one <- sum((seasonal_adjacencies$cold$RD==1)&(seasonal_adjacencies$cold$RL==1))
both_minus_one <- sum((seasonal_adjacencies$cold$RD==-1)&(seasonal_adjacencies$cold$RL==-1))
round(100*(both_one+both_minus_one)/N_links, digits=0)

# HOT HOT
# RA-RD 37
# RA-RL 34
# RD-RL 41

# HOT COLD
# RA-RD 23
# RA-RL 23
# RD-RL 26
# RA-RA 19
# RD-RD 26
# RL-RL 32

# COLD COLD
# RA-RD 37
# RA-RL 31
# RD-RL 35

#-------------------------------------------------------------------------------

#Components

get_graph_stats <- function(adj){
  
  g <- graph_from_adjacency_matrix(
    adj,
    mode = "undirected", weighted = TRUE, diag = FALSE, 
    add.colnames = NA, add.rownames = NULL
  )
  
  comp <- components(g)
  
  print(table(table(comp$membership)))
  
  #label_LCC <- names(table(comp$membership))[table(comp$membership)==max(table(comp$membership))]
  #nodes_LCC <- names(comp$membership)[comp$membership==label_LCC]
  #isolated_nodes <- names(table(comp$membership))[table(comp$membership)==1]
  
  #return(list("nodes_LCC"=nodes_LCC, "isolated_nodes"=isolated_nodes))
  
}

get_graph_stats(RA_hot_adj)

# RA cold
1   | 2  | 3   4   5  11 |  61 | 967 
771 | 55 | 8   4   3   1 |   1 |  1 

# RD cold
1   | 2  | 3    4    5    6 | 1092 
729 | 40 | 10   4    2    3 |   1 

# RL cold
1   | 2  | 3    4    5    8 | 1408 
471 | 29 | 7    1    1    1 |   1 

# RA hot
1   | 2  | 3    4    6 | 1362 
486 | 36 | 9    4    2 |  1 

# RD hot
1   | 2  | 3    4 | 1351 
511 | 35 | 13   1 |  1 

# RL hot
1   | 2  | 3    4    9 | 1366 
535 | 21 | 5    2    1 |  1 

#-------------------------------------------------------------------------------
