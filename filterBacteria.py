import pandas as pd
import numpy as np
from ete3 import NCBITaxa
ncbi = NCBITaxa()

# Questi li carico direttamente, ma sotto lascio il codice per fare vedere come ho ottenuto gli ids a partire dai nomi
# I nomi si ottengono come mostrato all'inizio dello script di R cleanData.r
species_names = pd.read_csv(path_species)
species_names = species_names[["names", "ncbi_ids"]]

#-------------------------------------------------------------------------------------------

# There are 109 missing ids, moreover some of them are WRONG, in the sense that the name and the ID don't correspond
#   Example: 
#       Line 4880 of species_names
#       Pseudomonas (a Bacteria) should have id 2052956 
#       Instead, its id here is 52956, that is actually the id of Spathodus erythrodon (a fish)
# pd.value_counts(species_names['ncbi_ids'].isnull())

# Therefore I throw the ids and to use the names for re-computing the ids via NCBI
# If you want, you can compare the recomputed ids and the original to see that some are wrong
ids_dict = ncbi.get_name_translator(species_names['names'])
ids_list = np.array(list(ids_dict.values()))[:,0]
name_list = np.array(ids_dict.keys())

# 192 names are NOT FOUND in ncbi database (maybe they're synonims to a more standard name)
# I can't use the corresponding id since it's not trustworthy
# Therefore, I just don't use them

# Questo è il file che ho caricato all'inizio
species_names = pd.DataFrame(
    {
        "name": name_list,
        "ncbi_id": ids_list
    }
)

#-------------------------------------------------------------------------------------------

superkingdom = []
for id in species_names['ncbi_id']:
    try:
        lineage = ncbi.get_lineage(id)
        superkingdom.append( ncbi.translate_to_names(lineage)[2] )
    except:
        superkingdom.append( 'not_found' )

species_names['superkingdom'] = superkingdom

#
# In superkingdom we get:
# Bacteria     4973
# Eukaryota    1329
# Archaea        48
#
# We only keep Bacteria
bacteria = species_names.loc[species_names['superkingdom']=='Bacteria']

# There are 5 ids that are repeated two times, which are: 
# 2585118, 208480, 103816, 332055, 121428
#
#                             name  ncbi_id
#354            Alistipes communis  2585118
#895          Bowdeniella nasicola   208480
#4685   Rhodococcus pyridinivorans   103816
#5090          Sphingobium indicum   332055
#5102       Sphingobium xenophagum   121428
#5931         Actinomyces nasicola   208480
#5936              Alistipes obesi  2585118
#6311  Rhodococcus biphenylivorans   103816
#6327    Sphingobium hydrophobicum   121428
#6328        Sphingobium japonicum   332055
#

# I will aggregate the counts in the following way:
#
# Actinomyces nasicola -> Bowdeniella nasicola
# Alistipes obesi -> Alistipes communis
# Rhodococcus biphenylivorans -> Rhodococcus pyridinivorans
# Sphingobium hydrophobicum -> Sphingobium xenophagum
# Sphingobium japonicum -> Sphingobium indicum
#

# Get ranks
id_to_rank = ncbi.get_rank(bacteria['ncbi_id'])
ranks = list(id_to_rank.values())
ids_with_rank = id_to_rank.keys()

# Keep only entries with rank information
bacteria = bacteria.loc[bacteria['ncbi_id'].isin(ids_with_rank), :]

bacteria['ncbi_id'] = ids_with_rank
bacteria['rank'] = ranks
# 
# species       4906
# strain          57
# subspecies       5
#

# Keep only species
bacteria = bacteria.loc[ bacteria['rank']=='species' ]

# Create aggregation table for subspecies and strains
# They will be aggregated into these species
# I add a column that is "aggregation_ids" that will be used later
aggregation_ids = np.array([])
rank_list = np.array([])
for i in bacteria.index:
    name = bacteria.loc[i, 'name']
    id = bacteria.loc[i, 'ncbi_id']
    rank = id_to_rank[id]
    rank_list = np.append(rank_list, rank)
    if rank!='species':
        lineage = ncbi.get_rank(ncbi.get_lineage(id))
        species_id = list(lineage.keys())[list(lineage.values()).index('species')]
        aggregation_ids = np.append(aggregation_ids, species_id)
    else:
        aggregation_ids = np.append(aggregation_ids, id)

aggregation_ids = aggregation_ids.astype(int)

# Non c'è qui il codice per la cosa dei sinonimi, l'ho fatto a mano sul terminale, ma ho aggregato anche loro

bacteria.loc[:, 'aggregation_id'] = aggregation_ids
bacteria.loc[:, 'rank'] = rank_list

# SAVE
bacteria = bacteria.reset_index()
bacteria = bacteria.drop('index', axis=1)
bacteria.to_csv("bacteria.csv", encoding='utf-8', index=True)
