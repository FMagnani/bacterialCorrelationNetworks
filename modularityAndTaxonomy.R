library(ggplot2)
library(data.table)
library(gridExtra)
library(ggpubr)

load(path_taxonomy_modularity)
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

translator <- function(taxlvl){
  
  if(taxlvl %in% c("phylum", "class", "order")){
    
    if(taxlvl=="phylum"){
      dict <- phylum_dict
    }else if(taxlvl=="class"){
      dict <- class_dict
    }else if(taxlvl=="order"){
      dict <- order_dict
    }
    
    id_to_name <- function(id){
      if(id %in% dict$ids){
        return(dict[dict$ids==id, ]$name)
      }else{
        return("Other")
      }
    }
    
  }else{
    id_to_name <- as.character
  }
  
  return(id_to_name)
  
}

#-------------------------------------------------------------------------------

max_modularity_RA <- sum(taxonomy_modularity$RA_adj$net_community$modularity)
max_modularity_RD <- sum(taxonomy_modularity$RD_adj$net_community$modularity)
max_modularity_RL <- sum(taxonomy_modularity$RL_adj$net_community$modularity)

max_single_modularity_RA <- max(taxonomy_modularity$RA_adj$net_community$modularity)
max_single_modularity_RD <- max(taxonomy_modularity$RD_adj$net_community$modularity)
max_single_modularity_RL <- max(taxonomy_modularity$RL_adj$net_community$modularity)

for(site in c("RA_adj", "RD_adj", "RL_adj")){
  for(lvl in names(taxonomy_modularity[[site]])){
    message(site, "    ", lvl, "    ", round(sum(taxonomy_modularity[[site]][[lvl]]$modularity), digits=3))
  }
}

#-------------------------------------------------------------------------------
  
get_taxlvl_df <- function(taxlvl){
  
  taxlvl_df <- rbind(
    data.frame(
      "tax"=taxonomy_modularity$RA[[taxlvl]]$tax,
      "modularity"=taxonomy_modularity$RA[[taxlvl]]$modularity,
      "site"="RA"
    ),
    data.frame(
      "tax"=taxonomy_modularity$RD[[taxlvl]]$tax,
      "modularity"=taxonomy_modularity$RD[[taxlvl]]$modularity,
      "site"="RD"
    ),
    data.frame(
      "tax"=taxonomy_modularity$RL[[taxlvl]]$tax,
      "modularity"=taxonomy_modularity$RL[[taxlvl]]$modularity,
      "site"="RL"
    )
  )
  
  return(taxlvl_df)
  
}

get_taxlvl_figure <- function(taxonomic_lvl, title){
  
  taxlvl_df <- get_taxlvl_df(taxonomic_lvl)
  taxlvl_df$tax_name <- unlist(lapply(taxlvl_df$tax, translator(taxonomic_lvl)))
  
  tax_to_show <- unique(taxlvl_df[taxlvl_df$modularity>0.005, ]$tax)
  taxname_to_show <- unlist(lapply(tax_to_show, translator(taxonomic_lvl)))
  
  mod_plt <- ggplot(
    taxlvl_df[taxlvl_df$tax %in% tax_to_show & taxlvl_df$tax_name!="-1" & taxlvl_df$tax_name!="Other", ]
  ) + geom_col(
    aes(
      x=modularity, 
      y=factor(tax_name, levels=unique(taxname_to_show)),
      fill=site
    ), position="dodge", color="black"
  ) + labs(
    title=title,
    x="", y=''
  ) + theme_minimal(
  )  + theme(
    axis.ticks.y=element_blank(), 
  )
  
  return(mod_plt)
  
}  

#-------------------------------------------------------------------------------

get_taxlvl_figure("phylum", "Phylum")

grid.arrange(
  top="Network modularity at various taxonomic levels",
  get_taxlvl_figure("phylum", "Phylum"),
  get_taxlvl_figure("class", "Class"),
  get_taxlvl_figure("order", "Order"),
  ncol=1
)
