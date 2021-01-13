library(vegan)
metadata=read.csv('~/Box Sync/EODF_analysis/raw_data/metadata/Metadata_16s_two_group.csv',row.names = 1)
asv=read.csv('~/Box Sync/EODF_analysis/intermediate_files/16s/dada2_output/GF only/asv_grouped.csv',header = T,row.names = 1)
samples_GF = rownames(metadata[metadata$Housing=='GF',])
samples_GF = gsub('-','.',samples_GF)
rownames(metadata)=samples_GF
samples_GF = intersect(samples_GF,colnames(asv))
asv= asv[,samples_GF]
metadata_GF = metadata[samples_GF,]
asv_relab = t(na.omit(data.frame(sapply(asv,function(x) x/sum(x)),row.names = rownames(asv))))
asv = t(asv)
asv.bc.dist=vegdist(asv_relab, method = "bray")
asv.bc.clust = hclust(asv.bc.dist, method = "average")

asv.bc.beta = betadisper(asv.bc.dist, metadata_GF$Group)
stat=adonis(asv.bc.dist~metadata_GF$Group, method = "bray")
scores(asv.bc.beta)
par(cex.lab=1.5)
par(cex.axis=1)
pcoa.all.fig = ordiplot(asv.bc.beta, type = "none", ylim = c(-0.4,0.4),xlim=c(-0.5,0.9),
                        xlab = paste("PCoA 1 (", round(asv.bc.beta$eig[1]/sum(asv.bc.beta$eig)*100,2), "% explained)"), 
                        ylab = paste("PCoA 2 (", round(asv.bc.beta$eig[2]/sum(asv.bc.beta$eig)*100,2), "% explained)"))
points(pcoa.all.fig, "sites", pch = 19, col = "blue", select = metadata_GF$Group == "WKY_Con_GF")
points(pcoa.all.fig, "sites", pch = 19, col = "purple", select = metadata_GF$Group == "WKY_EODF_GF")
points(pcoa.all.fig, "sites", pch = 19, col = "red", select = metadata_GF$Group == "SHRSP_Con_GF")
points(pcoa.all.fig, "sites", pch = 19, col = "green", select = metadata_GF$Group == "SHRSP_EODF_GF")
#ordiellipse(asv.bc.beta, metadata_GF$Group, conf = 0.95, label = F,col = c("red","green","blue","purple"))
ordispider(asv.bc.beta,metadata_GF$Group,col = c("red","green","blue","purple"))
legend("bottomright",legend = levels(metadata_GF$Group), col = c("red","green","blue","purple"), pch = 16)
mtext(paste("p = ", as.character(stat$aov.tab$`Pr(>F)`[1])))
dev.off()
rm(stat)
WKY_sample = rownames(metadata_GF[metadata_GF$Genotype == "WKY",])
SHRSP_sample = rownames(metadata_GF[metadata_GF$Genotype == "SHRSP",])
asv.wky.bc.dist = vegdist(asv_relab[WKY_sample,], method = "bray")
asv.wky.bc.clust = hclust(asv.wky.bc.dist, method = "average")
asv.shrsp.bc.dist = vegdist(asv_relab[SHRSP_sample,], method = "bray")
asv.shrsp.bc.clust = hclust(asv.shrsp.bc.dist, method = "average")
#stat.wky = anosim(asv.wky.bc.dist, grouping = metadata_GF[WKY_sample,"Group"], distance = "bray")
#stat.shrsp = anosim(asv.shrsp.bc.dist, grouping = metadata_GF[SHRSP_sample,"Group"], distance = "bray")
stat.wky=adonis(asv.wky.bc.dist~metadata_GF[WKY_sample,"Treatment"], method = "bray")
stat.shrsp=adonis(asv.shrsp.bc.dist~metadata_GF[SHRSP_sample,"Treatment"], method = "bray")
asv.wky.beta = betadisper(asv.wky.bc.dist, metadata_GF[WKY_sample,"Treatment"])
asv.shrsp.beta = betadisper(asv.shrsp.bc.dist, metadata_GF[SHRSP_sample,"Treatment"])

pcoa.wky.fig = ordiplot(asv.wky.beta, type = "none", ylim = c(-0.2,0.3), xlim = c(-0.2,0.3),
                        xlab = paste("PCoA 1 (", round(asv.wky.beta$eig[1]/sum(asv.wky.beta$eig)*100,2), "% explained)"), 
                        ylab = paste("PCoA 2 (", round(asv.wky.beta$eig[2]/sum(asv.wky.beta$eig)*100,2), "% explained)"))
points(pcoa.wky.fig, "sites", pch = 19, col = "blue", select = metadata_GF[WKY_sample,"Treatment"] == "Con")
points(pcoa.wky.fig, "sites", pch = 19, col = "purple", select = metadata_GF[WKY_sample,"Treatment"] == "EODF")
#ordiellipse(asv.wky.beta, metadata_GF[WKY_sample,"Treatment"], conf = 0.95,col = c("blue","purple"))
ordispider(asv.wky.beta,metadata_GF[WKY_sample,"Treatment"], col= c("blue","purple"))
legend("bottomright",legend = levels(metadata_GF$Treatment), col = c("blue","purple"), pch = 16)
mtext(paste("p = ", as.character(stat.wky$aov.tab$`Pr(>F)`[1])))
dev.off()

pcoa.shrsp.fig = ordiplot(asv.shrsp.beta, type = "none", ylim = c(-0.4,0.3), xlim = c(-0.2,0.7),
                          xlab = paste("PCoA 1 (", round(asv.shrsp.beta$eig[1]/sum(asv.shrsp.beta$eig)*100,2), "% explained)"), 
                          ylab = paste("PCoA 2 (", round(asv.shrsp.beta$eig[2]/sum(asv.shrsp.beta$eig)*100,2), "% explained)"))
points(pcoa.shrsp.fig, "sites", pch = 19, col = "red", select = metadata_GF[SHRSP_sample,"Treatment"] == "Con")
points(pcoa.shrsp.fig, "sites", pch = 19, col = "green", select = metadata_GF[SHRSP_sample,"Treatment"] == "EODF")
#ordiellipse(asv.shrsp.beta, metadata_GF[SHRSP_sample,"Treatment"], conf = 0.95,col = c("red","green"))
ordispider(asv.shrsp.beta,metadata_GF[SHRSP_sample,"Treatment"], col= c("red","green"))
legend("bottomright",legend = levels(metadata_GF$Treatment), col = c("red","green"), pch = 16)
mtext(paste("p = ", as.character(stat.shrsp$aov.tab$`Pr(>F)`[1])))
dev.off()
rm(stat.wky, stat.shrsp)

