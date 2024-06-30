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

# DATES
RL_df <- data.frame(
  "dates" = names(bacteria_species_counts$RL)[order(names(bacteria_species_counts$RL))],
  "frags" = colSums(bacteria_species_counts$RL)[order(names(bacteria_species_counts$RL))],
  "col" = "#619CFF", "site"="RL"
)
RA_df <- data.frame(
  "dates" = names(bacteria_species_counts$RA)[order(names(bacteria_species_counts$RA))],
  "frags" = colSums(bacteria_species_counts$RA)[order(names(bacteria_species_counts$RA))],
  "col" = "#F8766D", "site"="RA"
)
RD_df <- data.frame(
  "dates" = names(bacteria_species_counts$RD)[order(names(bacteria_species_counts$RD))],
  "frags" = colSums(bacteria_species_counts$RD)[order(names(bacteria_species_counts$RD))],
  "col" = "#00BA38", "site"="RD"
)

dates_df <- rbind(RA_df, RD_df, RL_df)

p_dates <- ggplot(
  dates_df
) + geom_col(
  aes(x = as.Date(dates, "%Y.%m.%d"), y = log(frags), fill = site),
  position="stack", color="black"
) + labs(
  title="Total number of fragment counts in samples",
  x="Date",
  y="Log of total fragment counts"
) + theme_minimal() + theme(
  text=element_text(size=20)
)
p_dates

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

# PARAMETER IDENTIFICATION - MEDIAN THRESHOLD

column <- "nonzero_medians"

RA_RD_median_cor <- cor(row_stats[row_stats$site=="RA", column], row_stats[row_stats$site=="RD", column])
RA_RL_median_cor <- cor(row_stats[row_stats$site=="RA", column], row_stats[row_stats$site=="RL", column])
RD_RL_median_cor <- cor(row_stats[row_stats$site=="RD", column], row_stats[row_stats$site=="RL", column])
message(RA_RD_median_cor)
message(RA_RL_median_cor)
message(RD_RL_median_cor)

avg_median_df <- data.frame("avg_median"=row_stats$nonzero_medians, "log_avg_median"=log(row_stats$nonzero_medians), "site"=row_stats$site)

median_groups_factor <- c("< 0", "0 - 5", "5 - 10", "10 - 100", "100 - 1000", "> 1000")

avg_median_df$Median <- "-"
avg_median_df$Median[avg_median_df$avg_median<=1] <- "< 0"
avg_median_df$Median[avg_median_df$avg_median>=1 & avg_median_df$avg_median<5] <- "0 - 5"
avg_median_df$Median[avg_median_df$avg_median>=5 & avg_median_df$avg_median<10] <- "5 - 10"
avg_median_df$Median[avg_median_df$avg_median>=10 & avg_median_df$avg_median<100] <- "10 - 100"
avg_median_df$Median[avg_median_df$avg_median>=100 & avg_median_df$avg_median<1000] <- "100 - 1000"
avg_median_df$Median[avg_median_df$avg_median>=1000] <- "> 1000"

avg_median_df$Median <- factor(avg_median_df$Median, levels=median_groups_factor)

p_RA <- ggplot(
  avg_median_df[avg_median_df$site=="RA" & avg_median_df$log_avg_median>0, ]
) + geom_histogram(
  aes(log_avg_median, fill=Median), color="black", binwidth=0.4
) + scale_fill_manual(
  values=c(
    "< 0" = "#A3123A",
    "0 - 5" = "#D83842",
    "5 - 10" = "#A1C893",
    "10 - 100" = "#63AC5E",
    "100 - 1000" = "#388248",
    "> 1000" = "#24693D"
  )
) + labs(
  title="RA",
  x="",
  y=""
) + theme_minimal() + theme(
  text=element_text(size=10)
)

p_RD <- ggplot(
  avg_median_df[avg_median_df$site=="RD" & avg_median_df$log_avg_median>0, ]
) + geom_histogram(
  aes(log_avg_median, fill=Median), color="black", binwidth=0.4
) + scale_fill_manual(
  values=c(
    "< 0" = "#A3123A",
    "0 - 5" = "#D83842",
    "5 - 10" = "#A1C893",
    "10 - 100" = "#63AC5E",
    "100 - 1000" = "#388248",
    "> 1000" = "#24693D"
  )
) + labs(
  title="RD",
  x="",
  y="Counts"
) + theme_minimal() + theme(
  text=element_text(size=10)
)

p_RL <- ggplot(
  avg_median_df[avg_median_df$site=="RL" & avg_median_df$log_avg_median>0, ]
) + geom_histogram(
  aes(log_avg_median, fill=Median), color="black", binwidth=0.4
) + scale_fill_manual(
  values=c(
    "< 0" = "#A3123A",
    "0 - 5" = "#D83842",
    "5 - 10" = "#A1C893",
    "10 - 100" = "#63AC5E",
    "100 - 1000" = "#388248",
    "> 1000" = "#24693D"
  )
) + labs(
  title="RL",
  x="Logarithm of the median",
  y=""
) + theme_minimal() + theme(
  text=element_text(size=10)
)

grid.arrange(
  top="Distribution of the typical numerosity of bacterial species",
  p_RA, p_RD, p_RL, ncol=1
)

#-------------------------------------------------------------------------------

# PARAMETER IDENTIFICATION - PREVALENCE THRESHOLD

p <- ggplot(
  row_stats
) + geom_histogram(
  aes(prevalence_and_numerosity, fill=site), position="dodge", color="black", binwidth=0.1
) + labs(
  title="Distribution of the frequency of the species",
  x="Percentage of samples in which the species has more than 10 fragment counts",
  y="Number of species"
) + scale_x_continuous(
  labels = scales::percent_format(accuracy = 1)
) + theme_minimal() + theme(
  text=element_text(size=20)
)
p

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

# table(table(c(RL_ids, RA_ids, RD_ids)))
# 1    2    3 
# 182  151 1975 

#-------------------------------------------------------------------------------

# VENN DIAGRAM RARE FILTERS
venn.diagram(
  x = list(
    ids[RA_prevnum_filter], ids[RA_median_filter]
  ), 
  category.names = c("Prevalence","Abundance"),
  filename="venn_rareFilter_RA.png",
  #print.mode="percent",
  output=FALSE, imagetype="tiff", resolution=600,
  col=c("#F8766D"),
  fill=c("navy", 'red4')
)

venn.diagram(
  x = list(
    ids[RD_prevnum_filter], ids[RD_median_filter]
  ), 
  category.names = c("Prevalence","Abundance"),
  filename="venn_rareFilter_RD.png",
  #print.mode="percent",
  output=FALSE, imagetype="tiff", resolution=600,
  col=c("#00BA38"),
  fill=c("navy", 'red4')
)

venn.diagram(
  x = list(
    ids[RL_prevnum_filter], ids[RL_median_filter]
  ), 
  category.names = c("Prevalence","Abundance"),
  filename="venn_rareFilter_RL.png",
  #print.mode="percent",
  output=FALSE, imagetype="tiff", resolution=600,
  col=c("#619CFF"),
  fill=c("navy", 'red4')
)

# Common species in sites
venn.diagram(
  x = list(
    RA_ids, RD_ids, RL_ids
  ), 
  category.names = c("RA","RD","RL"),
  filename="venn_percent.png",
  print.mode="percent",
  output=FALSE, imagetype="tiff", resolution=600,
  col=c("#F8766D", '#00BA38', '#619CFF'),
  fill=c("#F8766D", '#00BA38', '#619CFF')
)

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

# Zeros in the dataset

# Start:
#  RA 47.9 %
#  RD 49.2 %
#  RL 49.1 %

# Post rare filter:
# RA 8.1 %
# RD 9.1 %
# RL 7.9 %

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

# CORRELATIONS

RA_corr_test <- corr.test(
  x=t(RA_clr), 
  use="pairwise", method="pearson",
  adjust="bonferroni",
  ci=FALSE
)

RD_corr_test <- corr.test(
  x=t(RD_clr), 
  use="pairwise", method="pearson",
  adjust="bonferroni",
  ci=FALSE
)

RL_corr_test <- corr.test(
  x=t(RL_clr), 
  use="pairwise", method="pearson",
  adjust="bonferroni",
  ci=FALSE
)

RL_corr <- RL_corr_test$r
RA_corr <- RA_corr_test$r
RD_corr <- RD_corr_test$r
diag(RL_corr) <- 0
diag(RA_corr) <- 0
diag(RD_corr) <- 0

# Check that corr matrix are symmetric
#table(as.vector(
#  RL_corr[upper.tri(RL_corr)] == t(RL_corr)[upper.tri(t(RL_corr))]
#))

# Keep common species
RL_corr <- RL_corr[common_bacteria, common_bacteria]
RA_corr <- RA_corr[common_bacteria, common_bacteria]
RD_corr <- RD_corr[common_bacteria, common_bacteria]

correlations <- list(
  "RA" = RA_corr,
  "RD" = RD_corr,
  "RL" = RL_corr
)
save(correlations, file=path_correlations)

#-------------------------------------------------------------------------------

RA_pval <- RA_corr_test$p[common_bacteria, common_bacteria]
RD_pval <- RD_corr_test$p[common_bacteria, common_bacteria]
RL_pval <- RL_corr_test$p[common_bacteria, common_bacteria]

diag(RA_pval) <- 1
diag(RD_pval) <- 1
diag(RL_pval) <- 1

# This is because the adjusted p-values are on the upper triangular,
#  the not adjusted values on the lower triangular (worst design I've ever seen)
RA_pval[lower.tri(RA_pval)] <- t(RA_pval)[lower.tri(RA_pval)]
RD_pval[lower.tri(RD_pval)] <- t(RD_pval)[lower.tri(RD_pval)]
RL_pval[lower.tri(RL_pval)] <- t(RL_pval)[lower.tri(RL_pval)]

# Check symmetric
#table(as.vector(
#  RL_pval[upper.tri(RL_pval)] == t(RL_pval)[upper.tri(t(RL_pval))]
#))

#-------------------------------------------------------------------------------

# CORRELATION OF CORRELATIONS

corrcorr <- data.frame(
  "RL" = as.vector(as.matrix(RL_corr[upper.tri(RL_corr)])),
  "RA" = as.vector(as.matrix(RA_corr[upper.tri(RA_corr)])),
  "RD" = as.vector(as.matrix(RD_corr[upper.tri(RD_corr)]))
)

corr.test(corrcorr$RA, corrcorr$RD)
rl_ra_annotation <- "R: 0.71\np < 2.2e-16"
rl_rd_annotation <- "R: 0.75\np < 2.2e-16"
rd_ra_annotation <- "R: 0.73\np < 2.2e-16"

binwidth <- 0.01
linewidth <- 2

p_ra_rl <- ggplot(
  corrcorr, 
  aes(x=RA, y=RL)
) + geom_bin2d(
  binwidth=c(binwidth,binwidth)
) + scale_fill_gradient2(
  low = 'white', high='black',
  trans = "log10"
) + labs(
  title="R = 0.71",
  x="RA",
  y="RL"
) + geom_abline(  
  aes(slope = 1, intercept=0, col='red2'), 
  show.legend=TRUE, linewidth=linewidth
) + geom_abline(  
  aes(slope = 0.71, intercept=0, col='red4'),
  show.legend=TRUE, linewidth=linewidth, lty='dashed'
) + xlim(-1,1) + ylim(-1,1) + scale_colour_manual(
  name='',
  labels = c("Identity", "Linear Regression"), 
  values=c("red4", "red2")
) + theme(
  legend.position="bottom", text=element_text(size=15)
)

p_ra_rd <- ggplot(
  corrcorr, 
  aes(x=RA, y=RD)
) + geom_bin2d(
  binwidth=c(binwidth,binwidth)
) + scale_fill_gradient2(
  low = 'white', high='black',
  trans = "log10"
) + labs(
  title="R = 0.73",
  x="RA",
  y="RD"
) + geom_abline(  
  aes(slope = 1, intercept=0, col='red2'), 
  show.legend=FALSE, linewidth=linewidth
) + geom_abline(  
  aes(slope = 0.73, intercept=0, col='red4'),
  show.legend=FALSE, linewidth=linewidth, lty='dashed'
) + xlim(-1,1) + ylim(-1,1) + scale_colour_manual(
  name='',
  labels = c("Identity", "Linear Regression"), 
  values=c("red4", "red2")
) + theme(
  legend.position="bottom", text=element_text(size=15)
)

p_rd_rl <- ggplot(
  corrcorr, 
  aes(x=RD, y=RL)
) + geom_bin2d(
  binwidth=c(binwidth,binwidth)
) + scale_fill_gradient2(
  low = 'white', high='black',
  trans = "log10"
) + labs(
  title="R = 0.75",
  x="RD",
  y="RL"
) + geom_abline(  
  aes(slope = 1, intercept=0, col='red2'), 
  show.legend=FALSE, linewidth=linewidth
) + geom_abline(  
  aes(slope = 0.75, intercept=0, col='red4'),
  show.legend=FALSE, linewidth=linewidth, lty='dashed'
) + xlim(-1,1) + ylim(-1,1) + scale_colour_manual(
  name='',
  labels = c("Identity", "Linear Regression"), 
  values=c("red4", "red2")
) + theme(
  legend.position="bottom", text=element_text(size=15)
)

grid.arrange(
  top = text_grob(
    "Comparison of bacterial correlations in different sites",
    size = 18
  ),
  p_ra_rd, p_ra_rl, p_rd_rl, ncol=3
)

#cor(corrcorr$RD, corrcorr$RL)
# RA-RD: 0.73
# RL-RA: 0.71
# RL-RD: 0.75

#-------------------------------------------------------------------------------

# NETWORKS

density_to_correlation_threshold <- function(site_corr, site_pval, density_thr){
  
  N_links <- length(site_corr)
  index_thr <- 1 + as.integer(density_thr*N_links)
  
  significative_corr <- site_corr*apply(site_pval<0.05, c(1,2), as.integer)
  
  corr_thr <- abs(
    site_corr[
      order(abs(
        significative_corr
      ), decreasing=TRUE)
    ][[index_thr]]
  )
  
  unconnected_nodes <- 100*sum(rowSums((abs(significative_corr)>corr_thr))==0)/1975
  message(paste("unconnected nodes: ", unconnected_nodes, "%", sep=''))
  
  return(corr_thr)
  
}

RA_r_thr <- density_to_correlation_threshold(RA_corr, RA_pval, 1/100)
RD_r_thr <- density_to_correlation_threshold(RD_corr, RD_pval, 1/100)
RL_r_thr <- density_to_correlation_threshold(RL_corr, RL_pval, 1/100)

# MAKE ADJACENCIES

RA_adj <- RA_corr
RA_adj[(abs(RA_corr) < RA_r_thr)] <- 0
RA_adj[RA_pval > 0.05] <- 0
RA_adj[RA_adj > 0] <-  1
RA_adj[RA_adj < 0] <- -1

RD_adj <- RD_corr
RD_adj[(abs(RD_corr) < RD_r_thr)] <- 0
RD_adj[RD_pval > 0.05] <- 0
RD_adj[RD_adj > 0] <-  1
RD_adj[RD_adj < 0] <- -1

RL_adj <- RL_corr
RL_adj[(abs(RL_corr) < RL_r_thr)] <- 0
RL_adj[RL_pval > 0.05] <- 0
RL_adj[RL_adj > 0] <-  1
RL_adj[RL_adj < 0] <- -1

adjacencies <- list(
  "RA_adj" = RA_adj,
  "RD_adj" = RD_adj,
  "RL_adj" = RL_adj
)
save(adjacencies, file=path_adjacencies)

