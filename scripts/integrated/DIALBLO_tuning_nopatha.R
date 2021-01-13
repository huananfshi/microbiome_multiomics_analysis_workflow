source("~/Box Sync/EODF_analysis/scripts/utility/R_packages.R")
source("~/Box Sync/EODF_analysis/scripts/utility/read_metadata.R")
taxonomics = read.csv("~/Box Sync/EODF_analysis/intermediate_files/humann3/taxa_renorm.csv", header = TRUE,row.names = 1)
gene = read.csv("~/Box Sync/EODF_analysis/intermediate_files/humann3/ecs_cpm_all_nonspecies_name_renorm.csv",header = T, skip = 1,row.names = 1,numerals = 'no.loss',stringsAsFactors = F)
#pathway = read.csv("~/Box Sync/EODF_analysis/intermediate_files/humann3/pathways_cpm_all_nonspecies.csv",header = T,row.names = 1)
plasma = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/plasma_metabolites_by_subpwy.csv",header = T,row.names = 1)
cecal = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/cecal_metabolites_by_subpwy.csv",header = T,row.names = 1)

plasma = t(plasma)
cecal = t(cecal)
taxonomics=data.frame(t(taxonomics))
colnames(plasma) =lapply(colnames(plasma), function(x) paste('Plasma_',x,sep ="" ))
colnames(cecal) =lapply(colnames(cecal), function(x) paste('Cecal_',x,sep ="" ))

index = intersect(rownames(metadata),rownames(taxonomics))
index = intersect(index,rownames(gene))

#index = intersect(index,rownames(pathway))
index = intersect(index,rownames(plasma))
index = intersect(index,rownames(cecal))

data = list(taxonomy = taxonomics[index,],
            gene_abun = gene[index,],
            #pathway_abun = pathway[index,],
            plasma_metabolomics = plasma[index,],
            cecal_metabolomics =cecal[index,])
lapply(data,dim)
metadata = metadata[index,]
summary(metadata$Genotype_Group)


design = matrix(0.1, ncol = length(data), nrow = length(data), 
                dimnames = list(names(data), names(data)))
diag(design) = 0
#tuning
set.seed(123)
test.diablo = block.splsda(data,metadata$Genotype_Group,ncomp = 5,design = design)

perf.diablo = perf(test.diablo,validation = "Mfold",folds =5, nrepeat = 5)
plot(perf.diablo)
ncomp = perf.diablo$choice.ncomp$WeightedVote["Overall.BER", "centroids.dist"]
test.keepX = list(taxomony = c(seq(5,30,5)),
                  gene_abun=c(seq(15,50,10)),
                  #pathway_abun = c(seq(15,30,4)),
                  plasma_metabolomics = c(seq(5,30,4)),
                  cecal_metabolomics = c(seq(5,30,4))

tune.eodf = tune.block.splsda(X = data, Y = metadata$Genotype_Group,
                              ncomp = 2,test.keepX = test.keepX,
                              design = design,
                              validation = "Mfold",folds = 5,nrepeat = 5,
                              cpus = 4,dist = "centroids.dist")
list.keepX = tune.eodf$choice.keepX
save(data,metadata,file = "~/Box Sync/EODF_analysis/intermediate_files/humann3/mixOmics_nopatha.RDATA")
save(tune.eodf,design,ncomp,file="~/Box Sync/EODF_analysis/intermediate_files/humann3/DIABLO_tuning_nopatha.RDATA")

