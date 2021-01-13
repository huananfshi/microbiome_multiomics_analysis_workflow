set.seed(123)
#revised alpha diversity stats
counts = cbind(metadata,specnumber(species_filt))
colnames(counts)[4]='Count'
pbad2way(Count~Genotype*Group,data=counts)

shan = cbind(metadata,species_filt.alpha)
pbad2way(species_filt.alpha~Genotype*Group,data=shan)
sim = cbind(metadata,species_filt.sim)
pbad2way(species_filt.sim~Genotype*Group,data=sim)

even = cbind(metadata,species_filt.evenness)
pbad2way(species_filt.evenness~Genotype*Group,data=even)
pairwise.wilcox.test(counts$Count,counts$Genotype_Group,p.adjust.method = 'BH')
pairwise.wilcox.test(shan$species_filt.alpha,shan$Genotype_Group,p.adjust.method = 'BH')
pairwise.wilcox.test(sim$species_filt.sim,sim$Genotype_Group,p.adjust.method = 'BH')
pairwise.wilcox.test(even$species_filt.evenness,even$Genotype_Group,p.adjust.method = 'BH')

#revised taxonomic relative abundance stats
phyla = read.csv("~/Desktop/phyla_group.csv",header = T,row.names = 1)
phyla = data.frame(t(phyla))
phyla = phyla[rownames(metadata),]
phyla = cbind(metadata,phyla)
rm(list=ls())


pbad2way(Actinobacteria~Genotype*Group,data=phyla)
pairwise.wilcox.test(phyla$Actinobacteria,phyla$Genotype_Group,p.adjust.method = 'BH')

pbad2way(Bacteroidetes~Genotype*Group,data=phyla)
pairwise.wilcox.test(phyla$Bacteroidetes,phyla$Genotype_Group,p.adjust.method = 'BH')

pbad2way(Deferribacteres~Genotype*Group,data=phyla)
pairwise.wilcox.test(phyla$Deferribacteres,phyla$Genotype_Group,p.adjust.method = 'BH')

pbad2way(Firmicutes~Genotype*Group,data=phyla)
pairwise.wilcox.test(phyla$Firmicutes,phyla$Genotype_Group,p.adjust.method = 'BH')

pbad2way(Proteobacteria~Genotype*Group,data=phyla)
pairwise.wilcox.test(phyla$Proteobacteria,phyla$Genotype_Group,p.adjust.method = 'BH')

pbad2way(Verrucomicrobia~Genotype*Group,data=phyla)
pairwise.wilcox.test(phyla$Verrucomicrobia,phyla$Genotype_Group,p.adjust.method = 'BH')


para =cbind(metadata,species$`Parasutterella Parasutterella_excrementihominis`)
colnames(para)[4] = 'para'
shapiro.test(species$`Parasutterella Parasutterella_excrementihominis`)
pbad2way(para~Genotype*Group,data=para)
pairwise.wilcox.test(para$para,para$Genotype_Group,p.adjust.method = 'BH')

#revised qpcr abundance stats

pcrs = read.csv("~/Desktop/qpcrs.csv",header = T,row.names = 1)
pcrs = pcrs[rownames(metadata),]
pcrs = cbind(metadata,pcrs)

pbad2way(B_Il6~Genotype*Group,data=pcrs)
pairwise.wilcox.test(pcrs$B_Il6,pcrs$Genotype_Group,p.adjust.method = 'BH')

pbad2way(B_MCP1~Genotype*Group,data=pcrs)
pairwise.wilcox.test(pcrs$B_MCP1,pcrs$Genotype_Group,p.adjust.method = 'BH')

pbad2way(K_Il6~Genotype*Group,data=pcrs)
pairwise.wilcox.test(pcrs$K_Il6,pcrs$Genotype_Group,p.adjust.method = 'BH')

pbad2way(K_MCP1~Genotype*Group,data=pcrs)
pairwise.wilcox.test(pcrs$K_MCP1,pcrs$Genotype_Group,p.adjust.method = 'BH')


pcrs2 = read.csv("~/Desktop/qpcrs2.csv",header = T)


pbad2way(i_Tnfa~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$i_Tnfa,pcrs2$Genotype_Group,p.adjust.method = 'BH')
pbad2way(i_il17~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$i_il17,pcrs2$Genotype_Group,p.adjust.method = 'BH')
pbad2way(i_p47~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$i_p47,pcrs2$Genotype_Group,p.adjust.method = 'BH')

pbad2way(c_muc2~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$c_muc2,pcrs2$Genotype_Group,p.adjust.method = 'BH')



pbad2way(c_tlr2~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$c_tlr2,pcrs2$Genotype_Group,p.adjust.method = 'BH')
pbad2way(c_tlr4~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$c_tlr4,pcrs2$Genotype_Group,p.adjust.method = 'BH')

pbad2way(c_tnfa~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$c_tnfa,pcrs2$Genotype_Group,p.adjust.method = 'BH')

pbad2way(c_il1a~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$c_il1a,pcrs2$Genotype_Group,p.adjust.method = 'BH')

pbad2way(c_il17~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$c_il17,pcrs2$Genotype_Group,p.adjust.method = 'BH')

pbad2way(c_il1b~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$c_il1b,pcrs2$Genotype_Group,p.adjust.method = 'BH')

pbad2way(c_p47~Genotype*Group,data=pcrs2)
pairwise.wilcox.test(pcrs2$c_p47,pcrs2$Genotype_Group,p.adjust.method = 'BH')


#revised metagenomic stats
meta = read.csv("~/Desktop/metageno.csv",header = T)
pbad2way(X7adh~Genotype*Group,data=meta)
pairwise.wilcox.test(meta$X7adh,meta$Genotype_Group,p.adjust.method = 'BH')
pbad2way(bsh~Genotype*Group,data=meta)
pairwise.wilcox.test(meta$bsh,meta$Genotype_Group,p.adjust.method = 'BH')



pbad2way(bsh_bu~Genotype*Group,data=meta)
pairwise.wilcox.test(meta$bsh_bu,meta$Genotype_Group,p.adjust.method = 'BH')
pbad2way(bsh_ba~Genotype*Group,data=meta)
pairwise.wilcox.test(meta$bsh_ba,meta$Genotype_Group,p.adjust.method = 'BH')
pbad2way(bsh_lr~Genotype*Group,data=meta)
pairwise.wilcox.test(meta$bsh_lr,meta$Genotype_Group,p.adjust.method = 'BH')
pbad2way(bsh_lj~Genotype*Group,data=meta)
pairwise.wilcox.test(meta$bsh_lj,meta$Genotype_Group,p.adjust.method = 'BH')
