library(EnhancedVolcano)
plasma_res = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/plasma_metabolites_stat_for_R.csv",header = T,row.names = 1)
cecal_res = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/cecal_metabolites_stat_for_R.csv",header = T,row.names = 1)
#labeling BA
plasma_BA = read.csv("~/Box Sync/EODF_analysis/intermediate_files/BA analysis/plasma_ba.csv",header = T,row.names = 1)
cecal_BA = read.csv("~/Box Sync/EODF_analysis/intermediate_files/BA analysis/cecal_ba.csv",header = T,row.names = 1)

plasma_BA_label = rownames(plasma_BA)[1:21]
cecal_BA_label = rownames(cecal_BA)[1:17]

plasma_res[is.na(plasma_res)]=0
for(i in rownames(plasma_res)){
  if(plasma_res[i,"p_WC_SC"]==0.0000){
    plasma_res[i,"p_WC_SC"]=0.00001
  }
  
}


EnhancedVolcano(plasma_res,
                lab=rownames(plasma_res),
                x = "log2fc_SC_WC",
                y = "p_WC_SC",
                pCutoff = 10e-2,
                FCcutoff = 1,
                pLabellingCutoff = pCutoff,
                title="Plasma metabolites (SHRSP vs WKY)",
                pointSize=4,
                labSize=6,
                legendVisible=F,
                ylim=c(-0.5,5.5),
                selectLab = plasma_BA_label,
               border = "full",
               subtitle = "",
               drawConnectors = T,
               labhjust = 1,
               labvjust = 2,
               captionLabSize = 20,
               titleLabSize = 25
               
)


cecal_res[is.na(cecal_res)]=0
for(i in rownames(cecal_res)){
  if(cecal_res[i,"p_WC_SC"]==0.0000){
    cecal_res[i,"p_WC_SC"]=0.00001
  }
  
}
EnhancedVolcano(cecal_res,
                lab=rownames(cecal_res),
                x = "log2fc_SC_WC",
                y = "p_WC_SC",
                pCutoff = 10e-2,
                FCcutoff = 1,
                title="Cecal metabolites (SHRSP vs WKY)",
                pointSize=4,
                labSize=6,
                legendVisible=F,
                xlim = c(-6,6),
                border = "full",
                subtitle = "",
                selectLab = cecal_BA_label,
                ylim=c(-0.5,5.5),
                drawConnectors = T,
                labhjust = 1,
                labvjust = 2,
                captionLabSize = 20,
                titleLabSize = 25
)

#labeling BA
plasma_BA = read.csv("~/Box Sync/EODF_analysis/intermediate_files/BA analysis/plasma_ba.csv",header = T,row.names = 1)
cecal_BA = read.csv("~/Box Sync/EODF_analysis/intermediate_files/BA analysis/cecal_ba.csv",header = T,row.names = 1)

plasma_BA_label = rownames(plasma_BA)[1:21]
cecal_BA_label = rownames(cecal_BA)[1:17]

plasma_res[is.na(plasma_res)]=0
for(i in rownames(plasma_res)){
  if(plasma_res[i,"p_SC_SE"]==0.0000){
    plasma_res[i,"p_SC_SE"]=0.00001
  }
  
}


EnhancedVolcano(plasma_res,
                lab=rownames(plasma_res),
                x = "log2fc_SC_SE",
                y = "p_SC_SE",
                pCutoff = 10e-2,
                FCcutoff = 1,
                pLabellingCutoff = pCutoff,
                title="Plasma metabolites (SHRSP Con vs SHRSP EODF)",
                pointSize=4,
                labSize=6,
                legendVisible=F,
                xlim = c(-6,6),
                ylim=c(-0.5,5.5),
                selectLab = plasma_BA_label,
                border = "full",
                subtitle = "",
                drawConnectors = T,
                labhjust = 1,
                labvjust = 2,
                captionLabSize = 20,
                titleLabSize = 25
                
)


cecal_res[is.na(cecal_res)]=0
for(i in rownames(cecal_res)){
  if(cecal_res[i,"p_SC_SE"]==0.0000){
    cecal_res[i,"p_SC_SE"]=0.00001
  }
  
}
EnhancedVolcano(cecal_res,
                lab=rownames(cecal_res),
                x = "log2fc_SC_SE",
                y = "p_SC_SE",
                pCutoff = 10e-2,
                FCcutoff = 1,
                title="Cecal metabolites (SHRSP Con vs SHRSP EODF)",
                pointSize=4,
                labSize=6,
                legendVisible=F,
                xlim = c(-6,6),
                border = "full",
                subtitle = "",
                selectLab = cecal_BA_label,
                ylim=c(-0.5,5.5),
                drawConnectors = T,
                labhjust = 1,
                labvjust = 2,
                captionLabSize = 20,
                titleLabSize = 25
)

