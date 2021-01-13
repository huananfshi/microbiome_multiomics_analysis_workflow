source("~/Box Sync/EODF_analysis/scripts/utility/R_packages.R")
load("~/Box Sync/EODF_analysis/intermediate_files/humann3/mixOmics_nopatha.RDATA")
load("~/Box Sync/EODF_analysis/intermediate_files/humann3/DIABLO_tuning_nopatha.RDATA")
set.seed(123)
# list.keepX = list(gene_abun = c(13,50),
#                     pathway_abun = c(20,20),
#                     plasma_metabolomics = c(20,20),
#                     cecal_metabolomics = c(20,20)
# )
list.keepX = list(taxonomy = c(5,25),
                    gene_abun = c(15,35),
                    plasma_metabolomics = c(5,21),
                    cecal_metabolomics = c(29,29)
)
list.keepX = tune.eodf$choice.keepX

eodf.results = block.splsda(X=data,Y=metadata$Genotype_Group,ncomp = 2,
                            keepX=list.keepX, design = design)
eodf.perf = perf(eodf.results,dist=c("all", "max.dist", "centroids.dist", "mahalanobis.dist"),validation = "Mfold",folds = 5,nrepeat = 1, auc=F,progressBar = F,signif.threshold = 0.05)

eodf.test = DIABLO.test(eodf.results,method = "centroids.dist", validation = c("Mfold", "loo"), k = 5, nperm = 999, progress = TRUE)

eodf.features = list(taxonomy = selectVar(eodf.results, block = 'taxonomy', comp = 2)$taxonomy$name,
                     gene_abun = selectVar(eodf.results, block = 'gene_abun', comp = 2)$gene_abun$name,
                     pathway_abun = selectVar(eodf.results, block = 'pathway_abun', comp = 2)$pathway_abun$name,
                     plasma_metabolomics = selectVar(eodf.results, block = 'plasma_metabolomics', comp = 2)$plasma_metabolomics$name,
                     cecal_metabolomics = selectVar(eodf.results, block = 'cecal_metabolomics', comp = 2)$cecal_metabolomics$name 
) 




plotIndiv(eodf.results,ind.names = F,legend = T,
          subtitle = c("Taxonomy","Gene Abundance", "Pathway Abundance","Plasma Metabolomics","Cecal Metabolomics"),
          col = c("blue","purple","red","green"),
          pch=16,
          ellipse = T)
plotVar(eodf.results, var.names = F,legend = T,cutoff = 0.7)
par(cex.axis=1)
plotDiablo(eodf.results,col = c("blue","purple","red","green"),ncomp=1)
circosPlot(eodf.results, cutoff=0.7)
par(mar=c(1,1,1,1))
plotLoadings(eodf.results,block = "taxonomy",comp = 2,contrib = "max",title = "Taxonomy loading")
plotLoadings(eodf.results,block = "gene_abun",comp = 2,contrib = "max",title = "Gene abundance loading")
plotLoadings(eodf.results,block = "pathway_abun",comp = 2,contrib = "max", title = "Pathway abundance loading")
plotLoadings(eodf.results,block = "plasma_metabolomics",comp = 2,contrib = "max", title = "Plamsa metabolomics loading")
plotLoadings(eodf.results,block = "cecal_metabolomics",comp = 2,contrib = "max", title = "Cecal metabolomics loading")
# plotLoadings(eodf.results,block = "taxonomy",comp = 1,contrib = "max",title = "Taxonomy loading")
# plotLoadings(eodf.results,block = "pathway_abun",comp = 1,contrib = "max", title = "Pathway abundance loading")
# plotLoadings(eodf.results,block = "plasma_metabolomics",comp = 1,contrib = "max", title = "Plamsa metabolomics loading")
# plotLoadings(eodf.results,block = "cecal_metabolomics",comp = 1,contrib = "max", title = "Cecal metabolomics loading")
cimDiablo(eodf.results, legend.position = "right",transpose = T,margins = c(10,25))
# eodf.network = network(eodf.results,blocks = c(1,2,3,4),
#                        color.node = c("blue",'red', 'yellow', 'green'),
#                        cutoff = 0.6)
# write.graph(eodf.network$gR,file="~/Box Sync/EODF_analysis/outputs/DIABLO_network.gml",format= "gml")
auroc(eodf.results,roc.block = "taxonomy",roc.comp = 1, title = "Taxonomy ROC")
auroc(eodf.results,roc.block = "gene_abun",roc.comp = 1, title = "gene ROC")
auroc(eodf.results,roc.block = "pathway_abun",roc.comp = 1, title = "Pathway Abundance ROC" )
auroc(eodf.results,roc.block = "plasma_metabolomics",roc.comp = 1, title= "Plasma Metabolomics ROC")
auroc(eodf.results,roc.block = "cecal_metabolomics",roc.comp = 1,title= "Cecal Metabolomics ROC")
eodf.features$plasma_metabolomics = gsub("Plasma_","",eodf.features$plasma_metabolomics)
eodf.features$cecal_metabolomics = gsub("Cecal_","",eodf.features$cecal_metabolomics)
selectVar(eodf.results, block = 'gene_abun', comp = 1)$gene_abun$name
export(eodf.features,"~/Box Sync/EODF_analysis/outputs/humann3/mixOmics_features_selected_tuned_comp2.xlsx")
