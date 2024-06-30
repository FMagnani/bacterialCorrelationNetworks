library(data.table)
library(matrixStats)
library(readr)

#-------------------------------------------------------------------------------
# GET SPECIES NAMES AND IDS
# Non semplicissimo perché sono formattati maluccio
# I nomi sono ok, ma gli id per lo più sono sbagliati

txt <- read_file(path_abundances)
lines <- strsplit(txt, split='\n')[[1]]
row_names <- c()
#row_ids <- c()
for(l in lines){
  name <- strsplit(l, split=',')[[1]][[1]]
  #ncbi_id <- strsplit(strsplit(l, split=',')[[1]][[2]], split='\t')[[1]][[1]]
  row_names <- c(row_names, name)
  #row_ids <- c(row_ids, ncbi_id)
}

names <- row_names[-1] # First entry is header

#-------------------------------------------------------------------------------
# GET ABUNDANCES OF BACTERIA

abundances <- as.data.frame(fread(path_abundances))
rownames(abundances) <- names

# In questo file ottenuto con lo script di python filterBacteria.py ci sono solo
#   batteri e c'è un identificativo che indica come aggregare sinonimi, strains e subspecies
bacteria_names <- as.data.frame(fread(path_bacteria))

# Keep only bacteria
abundances <- abundances[bacteria_names$name, ]

#-------------------------------------------------------------------------------
# AGGREGATE COUNTS OF STRAINS, SUBSPECIES, SYNONIMS

abundances$aggregation_id <- bacteria_names$aggregation_id

aggregated_abundances <- aggregate(abundances[,-166], list(abundances$aggregation_id), FUN=sum)

new_rownames <- c()
for(aggr_id in aggregated_abundances$Group.1){
  n <- rownames(abundances[abundances$aggregation_id==aggr_id, ])[[1]]
  new_rownames <- c(new_rownames, n)
}

rownames(aggregated_abundances) <- new_rownames
aggregated_abundances$Group.1 <- NULL

abundances <- aggregated_abundances

#-------------------------------------------------------------------------------
# REMOVE ENTRIES WITHOUT METADATA
# I cannot know the site for them

is_entry_in_metadata <- function(entry_name){
  # Nello specifico mi riferisco alla colonna `Complete name` dei metadati
  # Questi nomi sono tipo "DTU_2020_1007150_1_MG_RA_CPH_Sewage_494"
  
  in_metadata <- TRUE
  
  # Se c'è un '-' non è nei metadati
  len <- length(strsplit(entry_name, split='-')[[1]])
  if(len!=1){
    in_metadata <- FALSE
  }

  # Se non inizia con 'DTU' non è nei metadati 
  init <- strsplit(entry_name, split='_')[[1]][[1]]
  if(init!="DTU"){
    in_metadata <- FALSE
  }
  
  return(in_metadata)
  
}

filter <- unlist(lapply(names(abundances), is_entry_in_metadata))
abundances <- abundances[ ,filter]

#-------------------------------------------------------------------------------
# RENAME ABUNDANCE COLUMNS

make_unique_id_metadata <- function(complete_name){
  # Make sample names like "DTU_2020_1007272_1_MG_RL"
  
  id <- strsplit(complete_name, split="_CPH")[[1]][[1]]
  return(id)
  
}

ids <- names(abundances)
names(abundances) <- unlist(lapply(ids, make_unique_id_metadata))

#-------------------------------------------------------------------------------
# SPLIT INTO SITE-SPECIFIC DATAFRAMES 

get_site_from_id <- function(id){
  return(
    tail(strsplit(id, split="_")[[1]], n=1)
  )
}

site_list <- unlist(lapply(names(abundances), get_site_from_id))
df_RA <- abundances[ ,site_list=="RA"]
df_RD <- abundances[ ,site_list=="RD"]
df_RL <- abundances[ ,site_list=="RL"]

#-------------------------------------------------------------------------------

metadata <- as.data.frame(fread(path_metadata))

# This is the same unique_id we gave to the columns of the abundance dataframe
metadata$'UNIQUE_ID' <- unlist(lapply(metadata$'Complete name', make_unique_id_metadata))

#-------------------------------------------------------------------------------
# RETRIEVE THE DATES FOR SITE-SPECIFIC DF 
# Rinomino le colonne tramite la data del sample corrispondente
# In questa fase filtro anche i replicati: se due sample hanno la stessa data
#   tengo quello con abbondanza maggiore

id_to_date <- function(id){
  date <- metadata$'Date of isolation'[metadata$UNIQUE_ID==id]
  return(date)
}

rename_site_specific_df <- function(df_site){
  # Put the date as colum names, also filtering the replicated measurements
  
  date_id_df <- data.frame(
    "date" = unlist(lapply(names(df_site), id_to_date)),
    "counts" = unlist(colSums(df_site))
  )
  
  names_to_keep <- c()
  # Di ogni data tengo il sample con più conteggi
  for(date in unique(date_id_df$date)){
    
    date_df <- date_id_df[date_id_df$date==date,]
    name_to_keep <- rownames(date_df)[date_df$counts==max(date_df$counts)]
    names_to_keep <- c(names_to_keep, name_to_keep)            
    
  }
  
  df_site <- df_site[,names_to_keep]
  names(df_site) <- unlist(lapply(names(df_site), id_to_date))
  
  return(df_site)
  
}

df_RA <- rename_site_specific_df(df_RA)
df_RD <- rename_site_specific_df(df_RD)
df_RL <- rename_site_specific_df(df_RL)

#-------------------------------------------------------------------------------
# SAVE

bacteria_species_counts <- list("RA"=df_RA, "RD"=df_RD, "RL"=df_RL)

save(bacteria_species_counts, file=path_bacteria_species_counts)
