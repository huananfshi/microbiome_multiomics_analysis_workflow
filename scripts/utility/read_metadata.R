#prerequest for all analysis
#read metadata, rows as samples, columns as features, csv format
metadata = read.csv("/Users/huananshi/Box Sync/EODF_analysis/raw_data/metadata/Bacteria Mapping File wo 30_32.csv", header = TRUE, row.names = 1)
rownames(metadata) = gsub(".BacteriaOnly_profile","",row.names(metadata))
rownames(metadata) = gsub("-",".",row.names(metadata))

head(metadata)
metadata$Genotype = ordered(metadata$Genotype, levels = c("WKY", "SHRSP"))
metadata$Group = ordered(metadata$Group, levels = c("Con", "EODF"))
metadata$Genotype_Group = ordered(metadata$Genotype_Group, levels = c("WKY_Con", "WKY_EODF","SHRSP_Con", "SHRSP_EODF"))
write.csv(metadata, 
          file = "~/Box Sync/EODF_analysis/intermediate_files/metadata/metadata.csv")
