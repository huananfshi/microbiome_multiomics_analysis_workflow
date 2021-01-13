#load libraries and prerequisite files 
source("~/Box Sync/EODF_analysis/scripts/utility/R_packages.R")
source("~/Box Sync/EODF_analysis/scripts/utility/read_metadata.R")
plasma = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/plasma_metabolites_log10.csv",header = T,row.names = 1,stringsAsFactors = F)
cecal = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/cecal_metabolites_log10.csv",header = T,row.names = 1,stringsAsFactors = F)
plasma_imp = read.csv("~/Box Sync/EODF_analysis/outputs/plasma_metabolites_RF_feature_importances_with_pwy_accuracy_100%_filtered.csv",header = T,stringsAsFactors = F)
cecal_imp = read.csv("~/Box Sync/EODF_analysis/outputs/cecal_metabolites_RF_feature_importances_with_pwy_accuracy_89%_filtered.csv",header = T,stringsAsFactors = F)
plasma_imp=plasma_imp[order(plasma_imp$importance,decreasing = T),]
cecal_imp=cecal_imp[order(cecal_imp$importance,decreasing = T),]
top_all_plasma = plasma_imp$BIOCHEMICAL[1:50]
top_all_cecal = cecal_imp$BIOCHEMICAL[1:50]
taxa = read.csv("~/Box Sync/EODF_analysis/intermediate_files/humann3/taxonomic_profiles_all_species.csv",header = T,row.names = 1)
taxa = data.frame(t(taxa))
species=rownames(taxa)
taxa = sapply(taxa, function(x) x/sum(x))
colSums(taxa)
taxa.renorm=asin(sqrt(taxa))
taxa.renorm=data.frame(t(taxa.renorm))
colnames(taxa.renorm)=species
taxa = taxa.renorm
rm(taxa.renorm)
colnames(taxa)=gsub('.*\\.','',colnames(taxa))

#write.csv(taxa,"~/Box Sync/EODF_analysis/intermediate_files/humann3/taxonomic_profiles_all_species_renorm.csv",quote = F)





plasma_index = colnames(plasma)
taxa_plasma = taxa[plasma_index,]

Df=metadata[colnames(plasma),]
plasma_top_imp = plasma_imp[1:50,]
plasma_top = plasma[rownames(plasma)%in%top_all_plasma,]
plasma_top_num = data.frame(lapply(plasma_top, function(x) as.numeric(x)),row.names = rownames(plasma_top))
plasma_top_num = data.frame(row.names = rownames(plasma_top))
corr.taxa.pla.met.top = cor(taxa_plasma,t(plasma_top_num),method = "spearman")
corr.taxa.pla.met.top = as.matrix(corr.taxa.pla.met.top)

col_anno = HeatmapAnnotation(Super.Pathway = plasma_top_imp$SUPER.PATHWAY,gp = gpar(brewer.pal(8,"Set3")), which = "column", importance = anno_barplot(plasma_top_imp$importance) )
#row_anno = rowAnnotation(df= Df$Genotype_Group, col = list(Group = c("SHRSP_Con" ="red","SHRSP_EODF" ="green","WKY_Con" = "blue","WKY_EODF" ="purple")),name = "Group", show_annotation_name = TRUE)

heatmap_plasma = Heatmap(corr.taxa.pla.met.top,
                         col=colorRampPalette(brewer.pal(8, "YlGnBu"))(25),
                         show_row_dend = F,
                         show_column_dend = FALSE,
                         show_column_names=F,
                         column_names_side = "bottom",
                         show_row_names=TRUE,
                         row_names_side = "left",
                         column_names_gp = gpar(fontsize = 10),
                         row_names_gp = gpar(fontsize = 10),
                         column_title="plasma metabolites",
                         column_title_side = "top",
                         column_title_gp = gpar(fontsize = 14),
                         clustering_distance_rows = "pearson",
                         clustering_method_rows = "average",
                         clustering_distance_columns = "pearson",
                         clustering_method_columns = "average",
                         heatmap_legend_param = list(at = c(-1,0,1),
                                                     labels = c("-1", "0", "1"),
                                                     color_bar = "continuous",
                                                     legend_direction="vertical",
                                                     labels_gp = gpar(fontsize = 8),
                                                     legend_height= unit(4, "cm"),
                                                     title_position = "leftcenter-rot",
                                                     title = "Spearman Correlation"
                         ),
                         show_heatmap_legend = TRUE,
                         top_annotation = col_anno,
                         height = unit(15, "cm"),
                         width = unit(10,'cm')
)
draw(heatmap_plasma,
     heatmap_legend_side = "right",
     row_title="Species",
     row_title_side="left",
     column_title = "",
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 8),
     row_title_gp = gpar(fontsize = 15)
)

cecal_index = colnames(cecal)
taxa_cecal = taxa[cecal_index,]
Df_cecal=metadata[colnames(cecal),]
cecal_top_imp = cecal_imp[1:50,]
cecal_top = cecal[rownames(cecal)%in%top_all_cecal,]

cecal_top_num = data.frame(lapply(cecal_top, function(x) as.numeric(x)) ,row.names = rownames(cecal_top))

corr.taxa.cec.met.top = cor(taxa_cecal,t(cecal_top_num),method = "spearman")
corr.taxa.cec.met.top = as.matrix(corr.taxa.cec.met.top)

col_anno_cecal = HeatmapAnnotation(Super.Pathway = cecal_top_imp$SUPER.PATHWAY,gp = gpar(brewer.pal(8,"Set3")), which = "column", importance = anno_barplot(cecal_top_imp$importance) )
heatmap_cecal = Heatmap(corr.taxa.cec.met.top,
                        col=colorRampPalette(brewer.pal(8, "YlGnBu"))(25),
                        show_row_dend = F,
                        show_column_dend = FALSE,
                        show_column_names=F,
                        column_names_side = "bottom",
                        show_row_names=TRUE,
                        row_names_side = "left",
                        column_names_gp = gpar(fontsize = 10),
                        row_names_gp = gpar(fontsize = 10),
                        column_title="cecal metabolites",
                        column_title_side = "top",
                        column_title_gp = gpar(fontsize = 14),
                        clustering_distance_rows = "pearson",
                        clustering_method_rows = "average",
                        clustering_distance_columns = "pearson",
                        clustering_method_columns = "average",
                        heatmap_legend_param = list(at = c(-1,0,1),
                                                    labels = c("-1", "0", "1"),
                                                    color_bar = "continuous",
                                                    legend_direction="vertical",
                                                    labels_gp = gpar(fontsize = 8),
                                                    legend_height= unit(4, "cm"),
                                                    title_position = "leftcenter-rot",
                                                    title = "Spearman Correlation"
                        ),
                        show_heatmap_legend = TRUE,
                        top_annotation = col_anno_cecal,
                        height = unit(15,'cm'),
                        width = unit(10, "cm"),
)

draw(heatmap_cecal,
     heatmap_legend_side = "right",
     row_title="Species",
     row_title_side="left",
     column_title = "",
     column_title_side = "bottom",
     column_title_gp = gpar(fontsize = 8),
     row_title_gp = gpar(fontsize = 15)
)

