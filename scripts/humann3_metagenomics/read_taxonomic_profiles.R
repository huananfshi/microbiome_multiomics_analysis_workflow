#read taxonomic_profiles
#format: rows as samples and columns as taxa, csv format, species level only
source("~/Box Sync/EODF_analysis/scripts/utility/R_packages.R")
source("~/Box Sync/EODF_analysis/scripts/utility/read_metadata.R")
taxa = read.csv(file = "~/Box Sync/EODF_analysis/raw_data/humann3/metaphlan/species.tsv", 
                header = TRUE, row.names = 1,  check.names = FALSE, sep='\t')

species = data.frame(taxa[2:ncol(taxa)])
head(row.names(species))

#clean row names to species only
rownames(species) = gsub(".*\\|g__", "", rownames(species))
rownames(species) = gsub("\\|s__", " ", rownames(species))

head(row.names(species))

#test sum = ~100
colSums(species)

#filter for beta diversity
dim(species[apply(species, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(species)), ])
species_filt = species[apply(species, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(species)), ]
head(row.names(species_filt))
Species_samples = intersect(rownames(metadata),colnames(species_filt))
metadata_taxa = metadata[Species_samples,] 
species_filt = data.frame(t(species_filt), check.names = F)
species_filt=species_filt[Species_samples,] 
species = data.frame(t(species), check.names = F)
species=species[Species_samples,] 
#define subgroups
#metadata is required before setup subgroups: for control only, index: metadata$Group == "Con" (default); 
#for SHRSP only index: metadata$Genotype == "SHRSP 
species_filt_con=species_filt[metadata_taxa$Group=="Con",]
species_con=species[metadata_taxa$Group=="Con",]
species_filt_shrsp=species_filt[metadata_taxa$Genotype=="SHRSP",]
species_shrsp=species[metadata_taxa$Genotype=="SHRSP",]
# #get genus level information for LEfSe
# # genus = taxa
# # tmp.ind = grep("\\|s__", rownames(genus), invert = T)
# # tmp = genus[tmp.ind,]
# 
# #Get species only information: this is an example for all four groups. 
# #Use different index to select subgroups
# #Use code in the comment section to delete information of certain levels
# # tmp.ind = grep("\\|g__", rownames(tmp))
# # genus = tmp[tmp.ind,]
# # rm(tmp,tmp.ind)
# # #rownames(genus) = gsub("k__.([^ ]*?\\|)", "", rownames(genus))
# # #rownames(genus) = gsub("c__.([^ ]*?\\|)", "", rownames(genus))
# # #rownames(genus) = gsub("o__.([^ ]*?\\|)", "", rownames(genus))
# # rownames(genus) = gsub(".*\\|", "", rownames(genus))
# # colnames(genus) = gsub("-",".",colnames(genus))
# # colnames(genus) = gsub("_taxonomic_profile","",colnames(genus))
# # genus = t(genus)
# # genus = genus[Species_samples,]
# # Genus_sample = as.data.frame(rownames(genus))
# # colnames(Genus_sample) = "Sample_ID"
# # #genus = cbind(Genus_sample,metadata_taxa$Genotype,metadata_taxa$Group,genus)
# # #colnames(genus)[2:3] = c("Genotype","Group")
# # genus = cbind(Genus_sample,metadata_taxa$Genotype_Group,genus)
# # colnames(genus)[2] = "Genotype_Group"
# # rownames(genus) = NULL
# 
# write.csv(species,
#           file = "~/Box Sync/EODF_analysis/intermediate_files/humann3/taxonomic_profiles_all_species.csv")
# species_con = cbind(metadata_taxa[metadata_taxa$Group=='Con','Genotype_Group'],species_con)
# colnames(species_con)[1]='Genotype_Group'
# write.table(species_con,sep = "\t", quote = F,
#           file ="~/Box Sync/EODF_analysis/intermediate_files/humann3/taxonomic_profiles_con_species.txt", row.names = FALSE)
# 
# species_shrsp = cbind(metadata_taxa[metadata_taxa$Genotype=='SHRSP','Genotype_Group'],species_shrsp)
# colnames(species_shrsp)[1]='Genotype_Group'
# write.table(species_shrsp,sep = "\t", quote = F,
#             file ="~/Box Sync/EODF_analysis/intermediate_files/humann3/taxonomic_profiles_shrsp_species.txt", row.names = FALSE)

#write.table(genus, sep = "\t", quote = F,
#            file = "~/Box Sync/EODF_analysis/intermediate_files/humann3/taxonomic_profiles_genus.txt", row.names = FALSE)
