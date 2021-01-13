source("~/Box Sync/EODF_analysis/scripts/utility/read_metadata.R")
source("~/Box Sync/EODF_analysis/scripts/utility/R_packages.R")
taxonomics = read.csv("~/Box Sync/EODF_analysis/intermediate_files/taxonomy/taxonomic_profiles_all_species.csv", header = TRUE,row.names = 1)
plasma = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/plasma_metabolites_by_subpwy.csv",header = T,row.names = 1)
cecal = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/cecal_metabolites_by_subpwy.csv",header = T,row.names = 1)
set.seed(123)
plasma = t(plasma)
cecal = t(cecal)

colnames(plasma) =lapply(colnames(plasma), function(x) paste('Plasma_',x,sep ="" ))
colnames(cecal) =lapply(colnames(cecal), function(x) paste('Cecal_',x,sep ="" ))
keepX_plasma = c(9,21)
keepX_cecal = c(9,9)
#Con
index_con = rownames(metadata[metadata$Group=="Con",])
index_con = intersect(index_con,rownames(plasma))
index_con = intersect(index_con,rownames(cecal))


metadata_con = metadata[index_con,]
summary(metadata_con$Genotype)

plasma_result_con = splsda(plasma[index_con,],metadata_con$Genotype,keepX = keepX_plasma)
cecal_result_con = splsda(cecal[index_con,],metadata_con$Genotype,keepX = keepX_cecal)
plotIndiv(plasma_result_con,ind.names = F,legend = F,ellipse = T,title = "Plasma Metabolomics",col = c("blue","red") ,pch = 16,X.label = "Component 1", Y.label = "Compoment 2",
          size.title = rel(3),size.xlabel = rel(2),size.ylabel = rel(2),size.axis = rel(1.5))
plotIndiv(cecal_result_con,ind.names = F,legend = F,ellipse = T,title = "Cecal Metabolomics",col = c("blue","red") ,pch = 16,X.label = "Component 1", Y.label = "Compoment 2",
          size.title = rel(3),size.xlabel = rel(2),size.ylabel = rel(2),size.axis = rel(1.5))

#SHRSP
index_shrsp = rownames(metadata[metadata$Genotype=="SHRSP",])
index_shrsp = intersect(index_shrsp,rownames(plasma))
index_shrsp = intersect(index_shrsp,rownames(cecal))


metadata_shrsp = metadata[index_shrsp,]
summary(metadata_shrsp$Group)

plasma_result_shrsp = splsda(plasma[index_shrsp,],metadata_shrsp$Group,keepX = keepX_plasma)
cecal_result_shrsp = splsda(cecal[index_shrsp,],metadata_shrsp$Group,keepX = keepX_cecal)
plotIndiv(plasma_result_shrsp,ind.names = F,legend = F,ellipse = T,title = "Plasma Metabolomics",col = c("red","green") ,pch = 16,X.label = "Component 1", Y.label = "Compoment 2",
          size.title = rel(3),size.xlabel = rel(2),size.ylabel = rel(2),size.axis = rel(1.5))
plotIndiv(cecal_result_shrsp,ind.names = F,legend = F,ellipse = T,title = "Cecal Metabolomics",col = c("red","green") ,pch = 16,X.label = "Component 1", Y.label = "Compoment 2",
          size.title = rel(3),size.xlabel = rel(2),size.ylabel = rel(2),size.axis = rel(1.5))

#four groups
index = rownames(metadata)
index = intersect(index,rownames(plasma))
index = intersect(index,rownames(cecal))


metadata = metadata[index,]
summary(metadata$Genotype_Group)

plasma_result = splsda(plasma[index,],metadata$Genotype_Group,keepX=c(9,21))
plasma.test =  MVA.test(plasma[index,],metadata$Genotype_Group,model = 'PLS-DA', ncomp=2,kout = 5, nperm = 999, progress = TRUE)
cecal_result = splsda(cecal[index,],metadata$Genotype_Group,ncomp = 2,keepX = c(9,9))
cecal.test =  MVA.test(cecal[index,],metadata$Genotype_Group,model = 'PLS-DA', ncomp=2,kout = 5, nperm = 999, progress = TRUE)

plotIndiv(plasma_result,ind.names = F,legend = T,ellipse = T,title = "Plasma Metabolomics",col = c("blue","purple","red","green") ,pch = 16,X.label = "Component 1", Y.label = "Compoment 2")
plotIndiv(cecal_result,ind.names = F,legend = T,ellipse = T,title = "Cecal Metabolomics",col = c("blue","purple","red","green") ,pch = 16,X.label = "Component 1", Y.label = "Compoment 2")
