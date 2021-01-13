#read gene family files (ECS)
source("~/Box Sync/EODF_analysis/scripts/utility/R_packages.R")
source("~/Box Sync/EODF_analysis/scripts/utility/read_metadata.R")
gene = read.csv(file = "~/Box Sync/EODF_analysis/intermediate_files/humann3/ecs_cpm_name.csv", header = TRUE, row.names = 1, check.names = FALSE)
head(rownames(gene)) 
#remove species information
tmp.ind = grep("\\|.*", rownames(gene), invert = T)
gene_unstratified = gene[tmp.ind,]
rm(tmp.ind) 

gene_name = sub("[^ ]*\\: {1}","",rownames(gene_unstratified))
rownames(gene_unstratified) = sub("\\:.*","", rownames(gene_unstratified))
head(rownames(gene_unstratified))


#test sum = ~1
colSums(gene_unstratified)


#filter for beta diversity 
dim(gene_unstratified[apply(gene_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(gene_unstratified)), ])
gene_unstratified_filt = gene_unstratified[apply(gene_unstratified, 1, function(x) sum(x > 0.0001) > 0.1 * ncol(gene_unstratified)), ]
Gene_samples = intersect(rownames(metadata),colnames(gene_unstratified_filt))
metadata_gene = metadata[Gene_samples,] 
gene_filt = data.frame(t(gene_unstratified_filt), check.names = F)
gene_filt = gene_filt[Gene_samples,]
gene = data.frame(t(gene_unstratified), check.names = F)
gene = gene[Gene_samples,]
rm(gene_unstratified,gene_unstratified_filt)

#define subgroups
#metadata is required before setup subgroups: for control only, index: metadata$Group == "Con" (default); 
#for SHRSP only index: metadata$Genotype == "SHRSP 
gene_filt_con=gene_filt[metadata$Group=="Con",]
gene_con=gene[metadata$Group=="Con",]

write.csv(gene,"~/Box Sync/EODF_analysis/intermediate_files/genes/ecs_relab_all_nonspecies.csv")
write.csv(gene_con,"~/Box Sync/EODF_analysis/intermediate_files/genes/ecs_relab_con_nonspecies.csv")

#add gene_names to ecs

gene_withname = rbind(gene_name,gene)
rownames(gene_withname)[1] = "Gene Name"
write.csv(gene_withname,"~/Box Sync/EODF_analysis/intermediate_files/humann3/ecs_cpm_all_nonspecies_name.csv")
