plasma_res = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/plasma_metabolites_stat_for_R.csv",header = T,row.names = 1)
cecal_res = read.csv("~/Box Sync/EODF_analysis/intermediate_files/metabolomics/cecal_metabolites_stat_for_R.csv",header = T,row.names = 1)

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
                title="Plasma metabolites (SHRSP vs WKY)",
                transcriptPointSize=4,
                transcriptLabSize=5,
                legendVisible=F,
                ylim=c(-0.5,5.5),
               border = "full",
               subtitle = ""
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
                transcriptPointSize=4,
                transcriptLabSize=5,
                legendVisible=F,
                xlim = c(-6,6),
                border = "full",
                subtitle = "",
                ylim=c(-0.5,5.5)
)
)
