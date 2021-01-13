#load libraries and prerequisite files 
source("~/Box Sync/EODF_analysis/scripts/utility/R_packages.R")
#source("~/Box Sync/EODF_analysis/scripts/utility/read_metadata.R")
source("~/Box Sync/EODF_analysis/scripts/humann3/read_taxonomic_profiles.R")
##Taxonomic analysis

#Alpha diversity
#richness/species count
#save plots to /outputs
png("~/Box Sync/EODF_analysis/outputs/humann3/Richness.png")
metadata_alpha = ordered(metadata_taxa$Genotype_Group,levels = c('WKY_Con','SHRSP_Con','WKY_EODF','SHRSP_EODF'))
boxplot(specnumber(species_filt)~metadata_alpha, ylab = "#number of species", xlab = "Genotype_Group", ylim = c(25,), outline = TRUE)
par(cex.lab=1.5)
par(cex.axis=1)
dev.off()
model.spec = lm(specnumber(species_filt)~metadata_taxa$Genotype + metadata_taxa$Group + metadata_taxa$Genotype:metadata_taxa$Group)
Anova(model.spec, type = "II" )
model.spec.aov = aov(specnumber(species_filt)~metadata_taxa$Genotype + metadata_taxa$Group + metadata_taxa$Genotype:metadata_taxa$Group)
TukeyHSD(model.spec.aov)
#kruskal.test(specnumber(species_filt[rownames(metadata_taxa),])~metadata_taxa$Genotype_Group)  
#shannon index
species_filt.alpha = vegan::diversity(species_filt, index = "shannon")
species_filt.sim = vegan::diversity(species_filt, index = "simpson")
species_filt.counts = apply(species_filt[,-1]>0,1,sum)
species_filt.evenness=species_filt.alpha/log(species_filt.counts)
png("~/Box Sync/EODF_analysis/outputs/humann3/Shannon.png")
boxplot(species_filt.alpha~metadata_alpha, xlab = "Genotype_Group", ylab = "Shannon Index", ylim = c(1,3))
dev.off()
model.alpha = lm(species_filt.alpha~metadata_taxa$Genotype+metadata_taxa$Group+metadata_taxa$Genotype:metadata_taxa$Group)
Anova(model.alpha, type = "II" )
model.alpha.aov = aov(species_filt.alpha~metadata_taxa$Genotype+metadata_taxa$Group+metadata_taxa$Genotype:metadata_taxa$Group)
TukeyHSD(model.alpha.aov)
#kruskal.test(species_filt.alpha~metadata_taxa$Genotype_Group)  
png("~/Box Sync/EODF_analysis/outputs/humann3/Simpson.png")
boxplot(species_filt.sim~metadata_alpha, xlab = "Genotype_Group", ylab = "Simpson Index", ylim = c(0.4,1))
dev.off()
model.alpha = lm(species_filt.sim~metadata_taxa$Genotype+metadata_taxa$Group+metadata_taxa$Genotype:metadata_taxa$Group)
Anova(model.alpha, type = "II" )
model.alpha.aov = aov(species_filt.sim~metadata_taxa$Genotype+metadata_taxa$Group+metadata_taxa$Genotype:metadata_taxa$Group)
TukeyHSD(model.alpha.aov)

png("~/Box Sync/EODF_analysis/outputs/humann3/Evenness.png")
boxplot(species_filt.evenness~metadata_alpha, xlab = "Genotype_Group", ylab = "Evenness", ylim = c(0.2,1))
dev.off()
model.alpha = lm(species_filt.evenness~metadata_taxa$Genotype+metadata_taxa$Group+metadata_taxa$Genotype:metadata_taxa$Group)
Anova(model.alpha, type = "II" )
model.alpha.aov = aov(species_filt.evenness~metadata_taxa$Genotype+metadata_taxa$Group+metadata_taxa$Genotype:metadata_taxa$Group)
TukeyHSD(model.alpha.aov)


#Beta Diversity
#distance calculation and clustering of taxonomy 
#bray-curtis
species_filt.bc.dist = vegdist(species_filt, method = "bray")
species_filt.bc.clust = hclust(species_filt.bc.dist, method = "average")
#dendrogram
#plot(species_filt.bc.clust)
#PCoA ordination
species_filt.beta = betadisper(species_filt.bc.dist, metadata_taxa$Genotype_Group)

#stat=anosim(species_filt.bc.dist, grouping = metadata_taxa$Genotype_Group, distance = "bray", permutations = 9999)
stat=adonis(species_filt.bc.dist~metadata_taxa$Genotype+metadata_taxa$Group + metadata_taxa$Genotype:metadata_taxa$Group, method = "bray")
#permutest(species_filt.beta,pairwise = FALSE, permutations = 9999)
scores(species_filt.beta)
#TukeyHSD(stat$)

#Plot PCoA ordination
png("~/Box Sync/EODF_analysis/outputs/humann3/PCoA_all.png")
par(cex.lab=1.5)
par(cex.axis=1)
pcoa.all.fig = ordiplot(species_filt.beta, type = "none", ylim = c(-0.6,0.6), 
                        xlab = paste("PCoA 1 (",round(species_filt.beta$eig[1]/sum(species_filt.beta$eig)*100,2), "% explained)"), 
                        ylab = paste("PCoA 2 (",round(species_filt.beta$eig[2]/sum(species_filt.beta$eig)*100,2), "% explained)"))
points(pcoa.all.fig, "sites", pch = 19, col = "blue", select = metadata_taxa$Genotype_Group == "WKY_Con")
points(pcoa.all.fig, "sites", pch = 19, col = "purple", select = metadata_taxa$Genotype_Group == "WKY_EODF")
points(pcoa.all.fig, "sites", pch = 19, col = "red", select = metadata_taxa$Genotype_Group == "SHRSP_Con")
points(pcoa.all.fig, "sites", pch = 19, col = "green", select = metadata_taxa$Genotype_Group == "SHRSP_EODF")
ordiellipse(species_filt.beta, metadata_taxa$Genotype_Group, conf = 0.95, label = F,col = c("blue","purple","red","green"))
#ordispider(species_filt.beta,metadata_taxa$Genotype_Group,col = c("blue","purple","red","green"))
legend("bottomright",legend = levels(metadata_taxa$Genotype_Group), col = c("blue","purple","red","green"), pch = 16)
mtext(paste("p = ", as.character(stat$aov.tab$`Pr(>F)`[1])))
dev.off()
# pcoa.all.fig.3d = ordiplot3d(species_filt.beta, type = "none")
# points(pcoa.all.fig.3d, "sites", pch = 19, col = "blue", select = metadata_taxa$Genotype_Group == "WKY_Con")
# points(pcoa.all.fig.3d, "sites", pch = 19, col = "purple", select = metadata_taxa$Genotype_Group == "WKY_EODF")
# points(pcoa.all.fig.3d, "sites", pch = 19, col = "red", select = metadata_taxa$Genotype_Group == "SHRSP_Con")
# points(pcoa.all.fig.3d, "sites", pch = 19, col = "green", select = metadata_taxa$Genotype_Group == "SHRSP_EODF")
rm(stat)

#Within genotype analysis 
WKY_sample = rownames(metadata_taxa[metadata_taxa$Genotype == "WKY",])
SHRSP_sample = rownames(metadata_taxa[metadata_taxa$Genotype == "SHRSP",])
species.wky.bc.dist = vegdist(species_filt[WKY_sample,], method = "bray")
species.wky.bc.clust = hclust(species.wky.bc.dist, method = "average")
species.shrsp.bc.dist = vegdist(species_filt[SHRSP_sample,], method = "bray")
species.shrsp.bc.clust = hclust(species.shrsp.bc.dist, method = "average")
#stat.wky = anosim(species.wky.bc.dist, grouping = metadata_taxa[WKY_sample,"Group"], distance = "bray")
#stat.shrsp = anosim(species.shrsp.bc.dist, grouping = metadata_taxa[SHRSP_sample,"Group"], distance = "bray")
stat.wky=adonis(species.wky.bc.dist~metadata_taxa[WKY_sample,"Group"], method = "bray")
stat.shrsp=adonis(species.shrsp.bc.dist~metadata_taxa[SHRSP_sample,"Group"], method = "bray")
species.wky.beta = betadisper(species.wky.bc.dist, metadata_taxa[WKY_sample,"Group"])
species.shrsp.beta = betadisper(species.shrsp.bc.dist, metadata_taxa[SHRSP_sample,"Group"])

png("~/Box Sync/EODF_analysis/outputs/humann/PCoA_wky.png")
pcoa.wky.fig = ordiplot(species.wky.beta, type = "none", ylim = c(-0.6,0.6), xlim = c(-0.5,0.5),
                        xlab = paste("PCoA 1 (", round(species.wky.beta$eig[1]/sum(species.wky.beta$eig)*100,2), "% explained)"), 
                        ylab = paste("PCoA 2 (", round(species.wky.beta$eig[2]/sum(species.wky.beta$eig)*100,2), "% explained)"))
points(pcoa.wky.fig, "sites", pch = 19, col = "blue", select = metadata_taxa[WKY_sample,"Group"] == "Con")
points(pcoa.wky.fig, "sites", pch = 19, col = "purple", select = metadata_taxa[WKY_sample,"Group"] == "EODF")
ordiellipse(species.wky.beta, metadata_taxa[WKY_sample,"Group"], conf = 0.95,col = c("blue","purple"))
#ordispider(species.wky.beta,metadata_taxa[WKY_sample,"Group"], col= c("blue","purple"))
legend("bottomright",legend = levels(metadata_taxa$Group), col = c("blue","purple"), pch = 16)
mtext(paste("p = ", as.character(stat.wky$aov.tab$`Pr(>F)`[1])))
dev.off()

png("~/Box Sync/EODF_analysis/outputs/humann3/PCoA_shrsp.png")
pcoa.shrsp.fig = ordiplot(species.shrsp.beta, type = "none", ylim = c(-0.6,0.6), xlim = c(-0.5,0.5),
                          xlab = paste("PCoA 1 (", round(species.shrsp.beta$eig[1]/sum(species.shrsp.beta$eig)*100,2), "% explained)"), 
                          ylab = paste("PCoA 2 (", round(species.shrsp.beta$eig[2]/sum(species.shrsp.beta$eig)*100,2), "% explained)"))
points(pcoa.shrsp.fig, "sites", pch = 19, col = "red", select = metadata_taxa[SHRSP_sample,"Group"] == "Con")
points(pcoa.shrsp.fig, "sites", pch = 19, col = "green", select = metadata_taxa[SHRSP_sample,"Group"] == "EODF")
ordiellipse(species.shrsp.beta, metadata_taxa[SHRSP_sample,"Group"], conf = 0.95,col = c("red","green"))
#ordispider(species.shrsp.beta,metadata_taxa[SHRSP_sample,"Group"], col= c("red","green"))
legend("bottomright",legend = levels(metadata_taxa$Group), col = c("red","green"), pch = 16)
mtext(paste("p = ", as.character(stat.shrsp$aov.tab$`Pr(>F)`[1])))
dev.off()
rm(stat.wky, stat.shrsp)


#Con analysis
Con_sample = rownames(metadata_taxa[metadata_taxa$Group=="Con",])
species.con.bc.dist = vegdist(species_filt_con, method = "bray")
species.con.bc.clust = hclust(species.con.bc.dist, method = "average")
#stat.con = anosim(species.con.bc.dist, grouping = metadata_taxa[Con_sample,"Genotype"], distance = "bray")
stat.con=adonis(species.con.bc.dist~metadata_taxa[Con_sample,"Genotype"], method = "bray")
species.con.beta = betadisper(species.con.bc.dist, metadata_taxa[Con_sample,"Genotype"])

png("~/Box Sync/EODF_analysis/outputs/humann3/PCoA_con.png")
pcoa.con.fig = ordiplot(species.con.beta, type = "none", ylim = c(-0.6,0.6), xlim = c(-0.5,0.5),
                        xlab = paste("PCoA 1 (", round(species.con.beta$eig[1]/sum(species.con.beta$eig)*100,2), "% explained)"), 
                        ylab = paste("PCoA 2 (", round(species.con.beta$eig[2]/sum(species.con.beta$eig)*100,2), "% explained)"))
points(pcoa.con.fig, "sites", pch = 19, col = "blue", select = metadata_taxa[Con_sample,"Genotype"] == "WKY")
points(pcoa.con.fig, "sites", pch = 19, col = "red", select = metadata_taxa[Con_sample,"Genotype"] == "SHRSP")
ordiellipse(species.con.beta, metadata_taxa[Con_sample,"Genotype"], conf = 0.95,col = c("blue","red"))
#ordispider(species.con.beta,metadata_taxa[Con_sample,"Genotype"], col= c("blue","red"))
legend("bottomright",legend = levels(metadata_taxa$Genotype), col = c("blue","red"), pch = 16)
mtext(paste("p = ", as.character(stat.con$aov.tab$`Pr(>F)`[1])))
dev.off()
rm(stat.con)

