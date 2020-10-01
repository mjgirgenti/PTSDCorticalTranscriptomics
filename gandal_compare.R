# Compare our DEG with gandall's
PLOT_PDF=F
areas <- area_by_sex
plots <- list()
plots_by_area <- list()
setwd("/SAY/standard2/IAC-CC0937-BiostatYSPH/dj333/.Project/PTSD/Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap")
for(ind in 1:length(area_by_sex)){
    ##L
    asd_meta = read.csv("./results/tables//Microarray_ASD_metaanalysis_092017.csv", row.names=1)
    scz_meta = read.csv("./results/tables/Microarray_SCZ_metaanalysis_092017.csv", row.names=1)
    bd_meta = read.csv("./results/tables/Microarray_BD_metaanalysis_092017.csv", row.names=1)
    mdd_meta = read.csv("./results/tables/Microarray_MDD_metaanalysis_092017.csv", row.names=1)
    aad_meta = read.csv("./results/tables/Microarray_AAD_metaanalysis_092017.csv", row.names=1)
    ibd_meta = read.csv("./results/tables/Microarray_IBD_metaanalysis_092017.csv", row.names=1)
    area <- areas[ind]
    load(deseq_files[ind])
    PTSD_meta <- sig02[, c('log2FoldChange', 'lfcSE', 'pvalue', 'padj', 'Genename', 'Geneid')]
    MDD_meta <- sig01[, c('log2FoldChange', 'lfcSE', 'pvalue', 'padj', 'Genename', 'Geneid')]    
    colnames(PTSD_meta) <- c('beta', 'SE', 'p', 'fdr', 'Genename', 'Geneid')
    colnames(MDD_meta) <- c('beta', 'SE', 'p', 'fdr', 'Genename', 'Geneid')    
    all_genes = intersect(intersect(intersect(intersect(intersect(intersect(rownames(asd_meta), rownames(scz_meta)), rownames(bd_meta)), rownames(aad_meta)), rownames(mdd_meta)), rownames(aad_meta)), PTSD_meta$Geneid)
    all_genes <- intersect(all_genes, MDD_meta$Geneid)
    
    allmeta = matrix(NA,nrow=length(all_genes), 8)
    allmeta[,1] = asd_meta$beta[match(all_genes, rownames(asd_meta))]
    allmeta[,2] = scz_meta$beta[match(all_genes, rownames(scz_meta))]
    allmeta[,3] = bd_meta$beta[match(all_genes, rownames(bd_meta))]
    allmeta[,4] = mdd_meta$beta[match(all_genes, rownames(mdd_meta))]
    allmeta[,5] = aad_meta$beta[match(all_genes, rownames(aad_meta))]
    allmeta[,6] = ibd_meta$beta[match(all_genes, rownames(ibd_meta))]
    allmeta[,7] = PTSD_meta$beta[match(all_genes, PTSD_meta$Geneid)]
    allmeta[,8] = MDD_meta$beta[match(all_genes, PTSD_meta$Geneid)]
    
    colnames(allmeta) = c("ASD", "SCZ", "BD", "MDD", "AAD", "IBD", paste0('PTSD/', area), paste0('MDD#1/', area))
    rownames(allmeta) = all_genes
    allmeta=  as.data.frame(allmeta)
    cor(allmeta,use="pairwise.complete.obs",method="spearman")
    
    null = data.frame(read.delim("./working_data//NullDistribution/null.txt",head=F))
    null = sort(null$V1)
    null = cbind(null, data.frame(prob=1-abs(2*seq(1:length(null))/length(null)-1)))
    
    
    #Make Bargraph
    comparisons = t(combn(seq(1,ncol(allmeta)),2))
    barplot = data.frame(Mean = NA, SEM=NA, p.fdr=NA)
    for (i in 1:dim(comparisons)[1]) {
      x = comparisons[i,1]
      y = comparisons[i,2]
      R = cor.test(allmeta[,x], allmeta[,y])
      rho =cor(allmeta[,x], allmeta[,y], method="spearman", use="pairwise.complete.obs")
      sem = (tanh(atanh(rho + 1.96/sqrt(nrow(allmeta)-3))) - tanh(atanh(rho - 1.96/sqrt(nrow(allmeta)-3))))/3.92
      
      barplot[i,] = c(rho, sem, R$p.value)
      rownames(barplot)[i] = paste(colnames(allmeta)[x],colnames(allmeta)[y],sep="-")
    }
    barplot$p <- barplot$p.fdr
    barplot$p.fdr = p.adjust(barplot$p.fdr,method="fdr")
    barplot$p.bootstrap = null$prob[findInterval(barplot$Mean, null$null)]
    barplot$p.symbol = ""
    barplot$p.symbol[barplot$p.bootstrap<0.05] = "*"
    barplot$p.symbol[barplot$p.bootstrap<0.01] = "**"
    barplot$p.symbol[barplot$p.bootstrap<0.001] = "***"
    barplot$Comparison = rownames(barplot)
    barplot$Sex <- sexcode1[[gsub("[0-9]", "", area_by_sex[ind])]]
    barplot$Area <- gsub("[A-Z]", "", area_by_sex[ind])
    barplot$modality="microarray"
    barplot <- barplot[grep('(PTSD|MDD#1)/', barplot$Comparison), ]
    barplot$comparison <- gsub('-(PTSD/[AFM][0-9]+|MDD#1/[AFM][0-9]+)', '', barplot$Comparison)
    num_plot <- length(plots_by_area[[gsub("[A-Z]", "", area_by_sex[ind])]])+1
    plots_by_area[[gsub("[A-Z]", "", area_by_sex[ind])]][[num_plot]] <- barplot
    Fig2_main = ggplot(barplot,aes(x = reorder(Comparison, -Mean), y=Mean, label=p.symbol)) + ylim(-1,1) +  
      geom_bar(stat="identity",fill="royalblue3",width=0.75) +
      geom_errorbar(aes(ymin=(Mean - SEM), ymax=(Mean + SEM)), position=position_dodge(width=0.8), width=0.25,size=0.25) +   
      #  ggtitle(title) +   
      theme(plot.title = element_text(size=20, face="bold", vjust=2)) +   
      labs(x="", y=expression(paste("Transcriptome correlation (", rho, ")", sep=""))) +     	
      theme(
        axis.text.x=element_text(angle=50, size=12, hjust=1),
        axis.text.y=element_text(size=10, vjust=0.5), 
        legend.title = element_text(size=12), 
        plot.title = element_text(size=12, face="bold"),
        axis.title.x = element_text(size=10, vjust=-0.35, face="bold"),
        axis.title.y = element_text(size=12, vjust=0.5),
        plot.margin=unit(c(2,2,1,2),"mm")
      ) + geom_text(color="red",size=4,aes(y=Mean+ sign(Mean)*SEM + sign(Mean)*.02))
    plots[[ind]] <- Fig2_main
    write.table(barplot, file = paste0(Figures_path, '/compare_with_gandall_bar_graph_', area_by_sex[ind], '.csv'), sep = ",", row.names = F)
}    
pdf(paste0(Figures_path, '/compare_with_gandall_bar_graph', '.pdf'))
plots
dev.off()

