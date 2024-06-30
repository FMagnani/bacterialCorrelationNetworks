library(ggplot2)
library(gridExtra)
library(ggalluvial)
library(stringr)

library(data.table)
library(matrixStats)
library(readr)
library(ggpubr)

taxtab <- as.data.frame(fread(path_taxonomy_table))

#-------------------------------------------------------------------------------

phylum_dict <- data.frame(
  "ids" = c(1224, 1239, 976, 201174, 200643),
  "names" = c('Pseudomonadota','Bacillota','Bacteroidota','Actinomycetota','Bacteroidia')
)

class_dict <- data.frame(
  "ids" = c(
    28211,1236,28216,186801,91061,1760,200643
  ),
  "names" = c(
    'Alphaproteobacteria', 'Gammaproteobacteria', 'Betaproteobacteria', 'Clostridia', 'Bacilli', 'Actinomycetes', 'Bacteroidia'
  )
)

order_dict <- data.frame(
  "ids" = c(186802,171549,80840,186826,91347,72274,2887326,204455,200644,204457,135614),
  "names" = c('Eubacteriales','Bacteroidales','Burkholderiales','Lactobacillales','Enterobacterales','Pseudomonadales','Moraxellales','Rhodobacterales','Flavobacteriales','Sphingomonadales','Xanthomonadales')
)

#-------------------------------------------------------------------------------

alluvial_data_long <- data.frame(
  "col1" = taxtab$phylum_id,
  "col2" = taxtab$class_id,
  "col3" = taxtab$order_id
)

alluvial_data_long$phylum_name <- "other"
for(i in  1:length(phylum_dict$ids)){
  alluvial_data_long[alluvial_data_long$col1==phylum_dict$ids[[i]], "phylum_name"] <- phylum_dict$names[[i]]
}
alluvial_data_long$class_name <- "other"
for(i in  1:length(class_dict$ids)){
  alluvial_data_long[alluvial_data_long$col2==class_dict$ids[[i]], "class_name"] <- class_dict$names[[i]]
}
alluvial_data_long$order_name <- "other"
for(i in  1:length(order_dict$ids)){
  alluvial_data_long[alluvial_data_long$col3==order_dict$ids[[i]], "order_name"] <- order_dict$names[[i]]
}

alluvial_data_long$flow <- paste(
  alluvial_data_long$phylum_name, 
  alluvial_data_long$class_name,
  alluvial_data_long$order_name,
  sep='_'
)
  
flow_table <- table(alluvial_data_long$flow)
alluvial_data_freq <- data.frame(
  "flow" = names(flow_table),
  "numerosity" = unname(flow_table)
)
alluvial_data_freq$numerosity.Var1 <- NULL
names(alluvial_data_freq) <- c("flow", "numerosity")
  
flow_df <- as.data.frame(str_split_fixed(alluvial_data_freq$flow, '_', 3))
names(flow_df) <- c("col1", "col2", "col3")
alluvial_data_freq[c("col1", "col2", "col3")] <- flow_df
  
p <- ggplot(
  data = alluvial_data_freq,
  aes(
    axis1=col1, axis2=col2, axis3=col3, y=numerosity
  )
) + scale_x_discrete(
  limits = c("Phylum", "Class", "Order"), 
  expand = c(.2, .05)
) + xlab("") + ylab("") + geom_alluvium(
  aes(fill=col1), width = 1/3.5, show.legend=FALSE
) + geom_stratum(
  width = 1/3.5, alpha=0.1, reverse = TRUE, color = "black", linewidth = 0.7
) + geom_text(
  stat = "stratum", aes(label = after_stat(stratum), fontface="bold"), size=4
) + theme_minimal() + ggtitle(
  "Taxonomic composition of the dataset"
) + theme(
  legend.title=element_text(size=15),
  axis.text.x = element_text(size=15), 
  axis.text.y = element_text(size=15),
  plot.title = element_text(size=18)
)
p  

