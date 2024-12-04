########################################################################################
############ UNM Herbarium ############
############ NM Checklist & Phylogeny 
############ Summarized: 20 Oct 2022
############ h.marx
############ Device: local
########################################################################################


require("taxize") #taxonomic resolution
require("ape") #read trees
require("dplyr")
require("rotl")


########################################## Current checklist of NM Flora: from Patrick Alexander
# read in data:
nmcheck <- read.csv(file = "data/NMChecklist/NMtax_13Oct22_UNM.csv")
nmcheck <- read_csv("NMtaxList_onSEINet.csv") #this is the checklist

head(nmcheck)
dim(nmcheck) #17057 species

length(unique(nmcheck$taxon)) #5809 taxa

# select genus + species: 
split <- strsplit(as.character(unique(nmcheck$taxon)), split=" ", fixed=TRUE)
ids.toremove <- sapply(split, function(i) length(i) < 2)
split <- split[!ids.toremove]
genus.name <- sapply(split, "[", 1L) #get just the genus_species_var...
species.name <- sapply(split, "[", 2L) #get just the species_var...
nmcheck_taxa <- unique(paste(genus.name, species.name, sep=" ")) #get  unique genus_species
length(unique(nmcheck_taxa)) #3972 taxa


########################################## Cleaned collections dataset (SEINet + GBIF): From Lizzie Lombardi
# read in data:
nmcol <- read.csv("MASTERDF_CLEANED.csv")
head(nmcol)
dim(nmcol) #258062

nmcol$scientificName

# select genus + species: 
split <- strsplit(as.character(nmcol$scientificName), split=" ", fixed=TRUE)
ids.toremove <- sapply(split, function(i) length(i) < 2)
split <- split[!ids.toremove]
genus.name <- sapply(split, "[", 1L) #get just the genus_species_var...
species.name <- sapply(split, "[", 2L) #get just the species_var...
nmcol_taxa <- unique(paste(genus.name, species.name, sep=" ")) #get unique genus_species
nmcol_taxa <- paste(genus.name, species.name, sep=" ")
length(nmcol_taxa) #8587 taxa


##########################################  Harmonize taxonomy across species lists (checklist + collections)

# Rotl: Match taxonomic names to the Open Tree Taxonomy.
resolved_names_nmcheck <- tnrs_match_names(nmcheck_taxa)
resolved_names_nmcol <- tnrs_match_names(nmcol_taxa)
dim(resolved_names_nmcheck) #3971; one duplicated
dim(resolved_names_nmcol) # 8514;  a lot of duplicated names

# add Rotl taxonomy resolution column #####!!!!!!! THIS IS CAUSING PROBLEMS
nmcheck_resolved <- add.taxized.column(df = nmcheck, colnum = 1, spliton = " ", sepas = "_") 
nmcol_resolved <- add.taxized.column(df = nmcol, colnum = 11, spliton = " ", sepas = "_")

dim(nmcol_resolved) #258062

# summarize counts of collections for resolved name ($taxized) 
head(nmcol_resolved)
nmcol_resolved %>% count(taxized, source, sort = T)
nmcol_resolved %>% count(taxized, sort = T) #harpo wants this list. Send. 
nmcol_resolved_count <- nmcol_resolved %>% count(taxized, sort = T)

# Combine:
resolved_names_nmcheck_unique <- unique(resolved_names_nmcheck$unique_name) # unique names in checklist
dim(resolved_names_nmcheck_unique) #3863
nmcol_resolved_count # number of collections per species from databases

resolved_names_nmcheck_unique <- as.data.frame(resolved_names_nmcheck_unique)
colnames(resolved_names_nmcheck_unique) <- "taxized_checklist"

nmdataset <- full_join(nmcol_resolved_count, resolved_names_nmcheck_unique, by = c("taxized" = "taxized_checklist"), keep=T)
head(nmdataset)
tail(nmdataset)
dim(nmdataset) #7499

dim(nmdataset %>% filter(!is.na(taxized))) #6817 databased collections from NM
dim(nmdataset %>% filter(is.na(taxized))) #682 checklist species not in databased collections 

dim(nmdataset %>% filter(!is.na(taxized_checklist)))#3862 NM checklist species
dim(nmdataset %>% filter(is.na(taxized_checklist))) #3637 collections from NM no in checklist 


########################################## Phylogeny(ies) from Smith & Brown 2018
# read in tree:
allotb <- read.tree("data/phylogeny/v0/ALLOTB.tre") #353185 tips
allotb$tip.label

nmcheck_taxa.phy <- gsub(nmdataset$taxized_checklist, pattern = " ", replacement = "_")
length(nmcheck_taxa.phy[(nmcheck_taxa.phy %in% allotb$tip.label)]) #3244
3244/3862 #0.8399793 taxa in phylogeny

nmcol_taxa.phy <- gsub(nmdataset$taxized, pattern = " ", replacement = "_")
length(nmcol_taxa.phy[(nmcol_taxa.phy %in% allotb$tip.label)]) #4137
4137/6817 #0.6068652 taxa in phylogeny

nmdataset$taxized[(!nmdataset$taxized %in% allotb$tip.label)] #species in collections that are not in phylogeny 
nmdataset$taxized_checklist[(!nmdataset$taxized_checklist %in% allotb$tip.label)] #species in checklist that are not in phylogeny 


########################################## Add columns for family and Order:





### TESTS ####
#########################################
# Could maybe use the OTL at some ponint...but for not keep it simple and just the Stephen's 
# https://lunasare.github.io/ssb2020_workshop/03-broken-taxa/index.html
resolve(nmcheck_taxa[1:10], db = "iplant")

iplant_resolve(sci=c("Helianthus annuus", "Homo sapiens"))

iplant_resolve("Helianthusss")

resolved_names_nmcheck <- tnrs_match_names(nmcheck_taxa)
head(resolved_names_nmcheck)
synonyms(resolved_names_nmcheck)
inspect(resolved_names_nmcheck, taxon_name = "Chenopodium cycloides")

resolved_names_nmcheck <- resolved_names_nmcheck %>% filter(!is.na(ott_id))

tol_about()
my_tree <- tol_induced_subtree(ott_ids = resolved_names_nmcheck$ott_id)
plot(my_tree)


ids.toremove <- sapply(split, function(i) length(i) < 2)

sapply(resolved_names_nmcheck, function(i) is_in_tree(resolved_names_nmcheck[i,]$ott_id))

is_in_tree(resolved_names[i,]$ott_id)

resolved_names_nmcheck[10:20,]

# taxonomic serial numbers (TSN) of the taxa from ITIS:
tsn <- get_tsn(nmcheck_taxa[10:20], accepted = FALSE)
# accepted names for each TSN:
lapply(tsn[-9], itis_acceptname) 


## https://github.com/iantrotter/ComTreeOpt/blob/master/MANUAL.md
devtools::install_github("iantrotter/ComTreeOpt")

library(ComTreeOpt)
library(picante)
?ComTreeOpt
data(FSN_species)
tree <- ComTreeOpt(FSN_species)
plot(tree, show.node.label=T)

data(EifelDataset_species)
tree <- ComTreeOpt(EifelDataset_species)
plot(tree, show.node.label=T)




