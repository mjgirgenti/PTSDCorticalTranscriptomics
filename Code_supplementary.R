#### codes for generating plots

### FIG 1a ====
sa <- paste0(rep(c("F","M","A"),each=4), rep(c(11,25,24,9),3))
files <- paste0("dat_ulval/deseq_",sa,"_m3_ulval_keep3.RData")
library(dplyr)
library(calibrate)
sa <- paste0(rep(c("F","M"),each=4), rep(c(11,25,24,9),2))
files <- paste0("dat_ulval/deseq_",sa,"_m3_ulval_keep3.RData")

pdf('dat_ulval/results/volcano_PTSD_deseq2_SA_0725.pdf', width=16, height=9)
par(mfrow=c(2,4), family="sans")
for (file in files){
  load(file)
  label <- strsplit(file,"_") %>% sapply("[",3) %>% gsub("r","",.)
  res <- sig02
  xmin <- min(res$log2FoldChange)-0.5
  xmax <- max(res$log2FoldChange)+0.5
  ymax <- max(-log10(res$pvalue)) + 0.5
  with(res, 
       plot(log2FoldChange, -log10(pvalue), pch=20, main=label, 
            xlim=c(xmin,xmax), ylim=c(0,ymax), xlab="log2(FoldChange)"))
  with(subset(res, padj<.05), 
       points(log2FoldChange, -log10(pvalue), pch=20, col="orangered"))
  #     abline(h=-log10(0.05), lty=2, col="green")
  if(min(res$padj)<.05)  
    abline(h=-log10(min(res$pvalue[res$padj>.05])), col="green", lty=2)
}
dev.off()
### FIG 1b ====
setwd('~/Documents/Projects/PTSD')
library(limma)
library(dplyr)
load(files[9]); sig11P <- sig02; sig11M <- sig01; a11 <- sig02$Geneid[sig02$padj<.05]
load(files[10]); sig25P <- sig02; sig25M <- sig01; a25 <- sig02$Geneid[sig02$padj<.05]
load(files[11]); sig24P <- sig02; sig24M <- sig01; a24 <- sig02$Geneid[sig02$padj<.05]
load(files[12]); sig9P <- sig02; sig9M <- sig01; a9 <- sig02$Geneid[sig02$padj<.05]
geneid <- union(a11,a25) %>% union(a24) %>% union(a9)

pdf('results/volcano&heatmaps/vennDiag_A_deseq2_0903.pdf')
sig11P <- sig11P[match(geneid,sig11P$Geneid),]
sig25P <- sig25P[match(geneid,sig25P$Geneid),]
sig24P <- sig24P[match(geneid,sig24P$Geneid),]
sig9P <- sig9P[match(geneid,sig9P$Geneid),]
com <- cbind(sig11P[c('Geneid','Genename','padj')],sig25P['padj'],sig24P['padj'],sig9P['padj'])
names(com) <- c("Geneid","Genename","A11","A25",'A24','A9')
com <- (com[,3:6]<.05)
for (k in 1:4) com[is.na(com[,k]),k] <- F
a <- vennCounts(com)
vennDiagram(com, circle.col = 1:4)
dev.off()

### FIG 1c ====
library(dplyr)
load('data/all_exon_junc.rdata')
load('data/degs_m3_ulval_keep3.rdata')
load('data/datProbes.RData')
exon_all <- datProbes[match(all_exon, datProbes$ensembl_gene_id),]
exon <- exon_all$external_gene_id[exon_all$chromosome_name %in% c(1:22,"X","Y") &
                                    exon_all$gene_biotype == "protein_coding"] %>% unique
all_junc[,1] <- all_junc[,1] %>% as.character
ct <- c()
for (i in 1:dim(all_junc)[1]){
  chr = all_junc$Chr[i]
  start = all_junc$Start[i]
  end = all_junc$End[i]
  id <- which(probes$chromosome_name==chr &
                probes$start_position<start &
                probes$end_position>end)
  if (length(id)) ct <- c(ct, id)
}
id_all <- unique(ct)
junc_all <- probes[id_all,]
junc <- junc_all$external_gene_name[junc_all$gene_biotype=="protein_coding" &
                                      junc_all$chromosome_name %in% c(1:22,"X","Y")]

##plot venn diagram
library(limma)
all_genes <- union(degs.A.05, exon) %>% union(., junc)
com <- data.frame(gene=all_genes %in% degs.A.05, 
                  exon=all_genes %in% exon,
                  junc=all_genes %in% junc)
a <- vennCounts(com)
pdf('results/pc&degs/gene_exon_junc.pdf')
vennDiagram(com, circle.col = 1:3)
dev.off()

### FIG 4c ====
library(limma)
pdf('results/pc&degs/vennDiag_FvM_deseq2_1126.pdf')
for (reg in c(11,25,24,9)){
  load(paste0('dat_ulval/deseq_F', reg, '_m3_ulval_keep3.RData'))
  sigF <- sig02
  load(paste0('dat_ulval/deseq_M', reg, '_m3_ulval_keep3.RData'))
  sigM <- sig02
  com <- merge(sigF, sigM, by="Geneid")[c(7,14)]< .05
  colnames(com) <- paste0(c("F","M"), reg)
  a <- vennCounts(com)
  vennDiagram(com, circle.col = 1:3)
}
dev.off()
pdf('results/pc&degs/vennDiag_MvP_deseq2_1126.pdf')
for (reg in c(11,25,24,9)){
  load(paste0('dat_ulval/deseq_A', reg, '_m3_ulval_keep3.RData'))
  com <- merge(sig01, sig02, by="Geneid")[c(7,14)]< .05
  colnames(com) <- paste0(c("M","P"), reg)
  a <- vennCounts(com)
  vennDiagram(com, circle.col = 1:3)
}
dev.off()

### FIG 4de ====
library(ggplot2)
degList <- list()
# cats <- 1:4 ##or use 12 including Axx
ncat <- 12
for (i in 1:ncat){
  load(files[i])
  genes <- as.character(sig02$Genename[sig02$padj<.05])
  genes <- unique(genes[genes!=""])
  degList[[sa[i]]] <- genes
  #   print(dim(sig02)[1])
}
ovlap <- matrix(ncol=ncat, nrow=ncat, NA)
for (i in 1:ncat) for (j in 1:ncat) ovlap[i,j] <- length(intersect(degList[[i]], degList[[j]]))
ovlapr <- ovlap
for (i in 1:ncat) for (j in 1:ncat) ovlapr[i,j] <- ovlap[i,j]/(ovlap[i,i]*ovlap[j,j])*15000
rownames(ovlapr)=sa[1:ncat]
colnames(ovlapr)=sa[1:ncat]
# heatmap(1-ovlapr, symm = T, Rowv = NA)
# cats <- c(5:12)  ##M
cats <- c(1:4,9:12) ##F
ovlap <- ovlap[cats,cats]
ovlapr <- ovlapr[cats,cats]
ovlapr[,3] <- 0; ovlapr[3,] <- 0; diag(ovlapr) <- 1 ##enrich
cormat <- ovlapr
cormat[lower.tri(cormat, diag=T)] <- NA
# Melt the correlation matrix
library(reshape2)
melted_cormat <- melt(cormat, na.rm = TRUE)

# Heatmap
library(ggplot2)
p1 <- ggplot(data = melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "lightyellow", 
                       midpoint = 0.5, limit = c(0,31), space = "Lab", 
                       name="Fold\nEnrichment") +
  #   theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()
# Get lower triangle of the correlation matrix
get_lower_tri<-function(cormat){
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}
# Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}
reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

# Reorder the correlation matrix
# cormat <- reorder_cormat(ovlapr)
cormat <- ovlapr
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)
# Create a ggheatmap
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "white", high = "red", mid = "lightyellow", 
                       midpoint = 0.5, limit = c(0,31), space = "Lab", 
                       name="Fold\nEnrichment") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


cntmat <- get_upper_tri(ovlap)
melted_cntmat <- melt(cntmat, na.rm=T)
pdf('dat_ulval/results/heatmap_overlap_FDR.05_F_0902.pdf')
ggheatmap + 
  geom_text(data = melted_cntmat, aes(Var2, Var1, label = value), color = "black", size = 4) +
  theme(
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.border = element_blank(),
    panel.background = element_blank(),
    axis.ticks = element_blank(),
    legend.justification = c(1, 0),
    legend.position = c(0.6, 0.8),
    legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,
                               title.position = "top", title.hjust = 0.5))
dev.off()


### FIG 5a ====
pdf('dat_ulval/results/vennDiag_PTSD_deseq2_SA_0703.pdf', width=10, height=5)

sigF <- sig11F; sigM <- sig11M
geneid <- intersect(sigF$Geneid, sigM$Geneid)
indF <- which(sigF$Geneid %in% geneid)
indM <- which(sigM$Geneid %in% geneid)
com <- left_join(sigF[indF,c('Geneid','padj')], sigM[indM,c('Geneid','padj')], by="Geneid")
names(com) <- c("Geneid","F11","M11")
com <- (com[,2:3]<.05)
df.vdc <- vennCounts(com)
class(df.vdc) <- "matrix"
df.vdc <- as.data.frame(df.vdc[c(3,2,4),])
df.vdc$x <- c(0,2,1)
df.vdc$y <- c(1,1,1)
df.venn <- data.frame(x = c(0, 2), y = c(1, 1),
                      labels = c("F11","M11"))
p1 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 0.5, color="white") +
  scale_fill_manual(values = c('orange', 'lightblue')) +
  coord_fixed() +
  theme_void() +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts, size = 5)

sigF <- sig25F; sigM <- sig25M
geneid <- intersect(sigF$Geneid, sigM$Geneid)
indF <- which(sigF$Geneid %in% geneid)
indM <- which(sigM$Geneid %in% geneid)
com <- left_join(sigF[indF,c('Geneid','padj')], sigM[indM,c('Geneid','padj')], by="Geneid")
names(com) <- c("Geneid","F25","M25")
com <- (com[,2:3]<.05)
df.vdc <- vennCounts(com)
class(df.vdc) <- "matrix"
df.vdc <- as.data.frame(df.vdc[c(3,2,4),])
df.vdc$x <- c(0,2,1)
df.vdc$y <- c(1,1,1)
df.venn <- data.frame(x = c(0, 2), y = c(1, 1),
                      labels = c("F25","M25"))
p2 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 0.5, color="white") +
  scale_fill_manual(values = c('orange', 'lightblue')) +
  coord_fixed() +
  theme_void() +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts, size = 5)

sigF <- sig24F; sigM <- sig24M
geneid <- intersect(sigF$Geneid, sigM$Geneid)
indF <- which(sigF$Geneid %in% geneid)
indM <- which(sigM$Geneid %in% geneid)
com <- left_join(sigF[indF,c('Geneid','padj')], sigM[indM,c('Geneid','padj')], by="Geneid")
names(com) <- c("Geneid","F24","M24")
com <- (com[,2:3]<.05)
df.vdc <- vennCounts(com)
class(df.vdc) <- "matrix"
df.vdc <- as.data.frame(df.vdc[c(3,2,4),])
df.vdc$x <- c(0,2,1)
df.vdc$y <- c(1,1,1)
df.venn <- data.frame(x = c(0, 2), y = c(1, 1),
                      labels = c("F24","M24"))
p3 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 0.5, color="white") +
  scale_fill_manual(values = c('orange', 'lightblue')) +
  coord_fixed() +
  theme_void() +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts, size = 5)

sigF <- sig9F; sigM <- sig9M
geneid <- intersect(sigF$Geneid, sigM$Geneid)
indF <- which(sigF$Geneid %in% geneid)
indM <- which(sigM$Geneid %in% geneid)
com <- left_join(sigF[indF,c('Geneid','padj')], sigM[indM,c('Geneid','padj')], by="Geneid")
names(com) <- c("Geneid","F9","M9")
com <- (com[,2:3]<.05)
df.vdc <- vennCounts(com)
class(df.vdc) <- "matrix"
df.vdc <- as.data.frame(df.vdc[c(3,2,4),])
df.vdc$x <- c(0,2,1)
df.vdc$y <- c(1,1,1)
df.venn <- data.frame(x = c(0, 2), y = c(1, 1),
                      labels = c("F9","M9"))
p4 <- ggplot(df.venn, aes(x0 = x, y0 = y, r = 1.5, fill = labels)) +
  geom_circle(alpha = .5, size = 0.5, color="white") +
  scale_fill_manual(values = c('orange', 'lightblue')) +
  coord_fixed() +
  theme_void() +
  labs(fill = NULL) +
  annotate("text", x = df.vdc$x, y = df.vdc$y, label = df.vdc$Counts, size = 5)

grid.arrange(p1,p2,p3,p4,nrow=2)
dev.off()

### FIG 5ce ====
genes <- xlsx::read.xlsx('data/GabaGenes forPlotting.xlsx', sheetIndex=1, header=F)
for (area in c(11,25,24,9)){
  load(paste0('dat_ulval/deseq_A', area, '_m3_ulval_keep3.RData'))
  temp1 <- sig01[match(genes$X1,sig01$Geneid),c('log2FoldChange','pvalue')]
  temp1$area <- area
  temp1$gene <- genes$X2
  temp2 <- sig02[match(genes$X1,sig02$Geneid),c('log2FoldChange','pvalue')]
  temp2$area <- area
  temp2$gene <- genes$X2
  if (area == 11){
    mdd <- temp1; ptsd <- temp2
  } else {
    mdd <- rbind(mdd, temp1); ptsd <- rbind(ptsd, temp2)
  }
}
mdd$area <- factor(mdd$area); ptsd$area <- factor(ptsd$area)
# library(gplots)
# balloonplot(mdd$gene, mdd$area, -log10(mdd$pvalue), dotsize= -log10(mdd$pvalue), label = F)
library(ggpubr)
mdd$`-logPvalue` = -log10(mdd$pvalue);
mdd$logPvalue <- log10(mdd$pvalue)
ptsd$`-logPvalue` = -log10(ptsd$pvalue);
ptsd$logPvalue <- log10(ptsd$pvalue)
pdf('results/pc&degs/GABA_logFC_pvalue.pdf', width=5.5, height=4)
ggballoonplot(mdd, x="gene", y="area", size="-logPvalue", fill="log2FoldChange") +
  gradient_fill(c("blue", "white", "red")) + labs(title="MDD GABA genes")
ggballoonplot(ptsd, x="gene", y="area", size="-logPvalue", fill="log2FoldChange") +
  gradient_fill(c("blue", "white", "red")) + labs(title="PTSD GABA genes")
dev.off()

##re-order
mdd$gene <- factor(mdd$gene, levels = c("ELFN1", "GAD2", "LHX6", "PNOC", "SLC32A1", 
                                        "SST", "GABRA1", "GABRA2", "GPHN", "PVALB", "VIP"))
ptsd$gene <- factor(ptsd$gene, levels = c("ELFN1", "GAD2", "LHX6", "PNOC", "SLC32A1", 
                                          "SST", "GABRA1", "GABRA2", "GPHN", "PVALB", "VIP"))
pdf('results/pc&degs/GABA_logFC_pvalue_1216.pdf', width=5.5, height=4)
ggballoonplot(mdd, x="gene", y="area", size="-logPvalue", fill="log2FoldChange") +
  gradient_fill(c("blue", "white", "red")) + labs(title="MDD GABA genes")
ggballoonplot(ptsd, x="gene", y="area", size="-logPvalue", fill="log2FoldChange") +
  gradient_fill(c("blue", "white", "red")) + labs(title="PTSD GABA genes")
dev.off()



### FIG S2b, FIG S7b ====
library(ggplot2)
library(xlsx)
rm(list=ls())
dt2 <- read.xlsx('results/wgcna&csea/WGCNA_ulval_qreg_CP_reg2_go.xlsx', sheetName="M")
dt2 <- dt2[order(dt2$module,dt2$p_value),]
dt2 <- dt2[dt2$source!="TF",]
dt2 <- dt2[dt2$module %in% names(table(dt2$module)[table(dt2$module)>1]),]
dt2$num <- sapply(1:nrow(dt2), function(n) sum(dt2$module[1:n]==dt2$module[n]))
dt2 <- dt2[dt2$num %in% c(1:2),]
dt2$minP <- sapply(dt2$module, function(m) min(dt2$p_value[dt2$module==m]))
dt2 <- dt2[order(dt2$minP, dt2$num),]
dt2$`-log10FDR` <- -log10(dt2$p_value)
dt2$id <- paste0("mod", as.numeric(dt2$module), "_", dt2$term_name)
dt2$id <- factor(dt2$id, levels=rev(as.character(dt2$id)))
dt2$term_name1 <- sapply(dt2$term_name, function(t) paste(strwrap(t,width=80),collapse="\n"))

pdf('results/wgcna&csea/WGCNA_CP_GO_M_replot.pdf', width=12, height=10)
ggplot(data=dt2, aes(x=id, y=`-log10FDR`, fill="black")) +
  geom_point(stat='identity', size=4, aes(col=module))  +
  scale_color_manual(values=levels(factor(dt2$module))) +
  # geom_segment(aes(y = 0, x = id, yend = `-log10FDR`, xend = id), color = "black") +
  xlab(NULL)+ylab("-log10(P_adj)")+geom_hline(yintercept = -log10(0.05), lty=2)+
  theme_minimal()+coord_flip()
dev.off()

### FIG S3 ====
rm(list=ls())
setwd('~/Documents/Projects/PTSD')
library(ggplot2)
library(dplyr)
library(gplots)
library(WGCNA)
labs = c("A","F","M"); labels <- c("Female+Male","Female","Male")
files <- paste0('dat_ulval/WGCNA_', labs, '_ulval_qreg_CP_reg2.RData')

pdf('results/wgcna&csea/module_csea_sig_0326.pdf', width=8, height=8)
for (i in 1:3){
  load(files[i])
  label = labels[i]
  to_plot <- to_plot[apply(to_plot,1,min)<.05,]
  to_plot <- to_plot[order(rownames(to_plot)),]
  
  heatmap.2(-log10(to_plot),
            col=blueWhiteRed(1000,1)[500:1000],
            scale="none",
            trace="none",
            cexRow = .5,
            cexCol = .5, 
            density.info = "none",
            colsep=0:6,
            rowsep=0:dim(to_plot)[1],
            sepcolor="grey",
            sepwidth=c(0.05,0.05),
            srtCol=45,
            offsetRow=0,
            offsetCol=-0.5,
            dendrogram="none",
            Rowv="none",
            Colv="none",
            key=T,
            key.xlab="-log10(P)",
            cellnote=signif(to_plot,3),
            notecex=.6,
            notecol="black",
            margins = c(5,8),
            main=paste0(label, " CSEA"))
}
dev.off()


### FIG S4 ====
res <- read.csv('results/cibersort/cibersortx/CIBERSORTx_Job9_Results.csv')
res1 <- res[c(1,5,12,22)]
res1$Astro <- rowSums(res[c(7,15,18,25)])
res1$ExN <- rowSums(res[c(2:4,6,8,14,19:21,24)])
res1$IntN <- rowSums(res[c(9,11,13,16:17,26:28,30)])
res1$OPC <- rowSums(res[c(10,29)])
resp <- melt(res1)
ncell <- 7
meds <- apply(res1[,1+1:ncell],2,median)
cells <- names(res1)[1+1:ncell]
resp$variable <- factor(resp$variable, levels=cells[order(meds, decreasing = T)])
pdf('~/Documents/Projects/PTSD/results/cibersort/cibersortx/fractions_cx_c29_comb2.pdf', width=12)
## comb2: S-mode
p <- ggplot(resp, aes(x=variable, y=value)) + 
  geom_boxplot() + 
  labs(title = "Cell type fractions based on snRNAseq") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
print(p)
dev.off()

library(reshape2)
library(ggplot2)
library(ggpubr)
library(RColorBrewer)
cmp <- melt(data = res1,id.vars = "Mixture") %>% as.data.frame
names(cmp) <- c("sampleID","ctype","proportion")
load('~/Documents/Projects/PTSD/data/dat13_ulval_keep3.RData')
ord <- match(as.character(cmp$sampleID),meta13uv_keep$Sample)
cmp$group <- meta13uv_keep$PrimaryDx[ord] %>% factor
cmp$region <- meta13uv_keep$Region[ord] %>% factor
cmp$sex <- meta13uv_keep$Sex[ord] %>% factor
cmp.ori <- cmp
cmp$ctype <- factor(cmp$ctype, levels=c("ExN","Astro","IntN","OPC","Endo","Oligo","Microglia"))
View(cmp)

ncell = 7
for (reg in c(9,11,24,25)){
  comp2 <- data.frame(matrix(NA, nrow=ncell, ncol=6))
  for(i in 1:ncell){
    c = levels(factor(cmp$ctype))[i]
    for(j in 1:3){
      sel1 <- cmp$group!=levels(cmp$group)[j] & cmp$region == reg & cmp$ctype==c
      res <- aov(proportion~group, data=cmp[sel1,]) %>% summary
      comp2[i,1:2+2*(j-1)] <- c(res[[1]]$`F value`[1], res[[1]]$`Pr(>F)`[1])
    }}
  rownames(comp2) <- levels(cmp$ctype)
  colnames(comp2) <- paste0(rep(c("MDDvsPTSD","CONvsPTSD","CONvsMDD"),each=2), "_", rep(c("F.value","P.value"),3))
  View(comp2)
  write.csv(comp2, file=paste0("results/cibersort/cibersortx/cibersortx_pvals_ba",reg,"_by_dx_comb2.csv"))
}


### FIG S5 ====
library(dplyr)
rm(list=ls())
load('data/dat_ulval_A_qnorm_noM.RData')
fpkm <- datExpr
demos <- datMeta
index=1:dim(datMeta)[1]#which(demos$Sex=="M")
seltpm=fpkm[,index]
seldem=datMeta[index,]
pca.tpm <- prcomp(t(seltpm))
pcatpm <- pca.tpm$x
datMeta$Sex = as.numeric(datMeta$Sex)
datMeta$PTSD = as.numeric(datMeta$PrimaryDx=="PTSD") %>% factor
datMeta$MDD = as.numeric(datMeta$PrimaryDx=="MDD") %>% factor
datMeta$Region = as.factor(datMeta$Region)
datMeta$Race[datMeta$Race=="AA"] <- "B"
datMeta$Race[datMeta$Race=="CAUC"] <- "W"
datMeta$Race <- factor(datMeta$Race)

features <- c('PrimaryDx','Sex','AgeDeath','RIN1','PMI','Region','Race','Center')
features_id <- match(features, names(datMeta))
npc <- 10
rf <- sapply(1:8, function(i){
  lm.pc <- lm(pcatpm[,1:npc]~datMeta[,features_id[i]])
  lapply(summary(lm.pc), function(x) x$r.squared) %>% unlist
})
colnames(rf) <- features 
rownames(rf) <- paste0('PC',1:npc)
write.csv(rf, file="results/pc&degs/pcs_rs_pc10.csv")

library("RColorBrewer")
library(gplots)
library(ggplot2)
pdf('results/pc&degs/rs_pc10_traits_replot.pdf')
colfunc <- colorRampPalette(c("blue", "white", "red"))
heatmap.2(as.matrix(rf), dendrogram = 'none', trace="none", Rowv = NULL, Colv=NULL, 
          density.info="none", col=brewer.pal(9,'YlOrRd'), cexCol = .9, cexRow = .9)
dev.off()

pvar <- pca.tpm$sdev; pvar <- pvar^2/sum(pvar^2)
prf <- colSums(pvar[1:npc]/sum(pvar[1:npc])*rf)
df.prf <- data.frame(trait=names(prf), var_explained=prf)
df.prf$trait <- factor(df.prf$trait, levels=df.prf$trait[order(df.prf$var_explained)])

## pie chart
pdf('results/pc&degs/rs_pc10_traits_replot3.pdf', width=10, height=8)
df.prf$traitPlus <- paste0(df.prf$trait, ": ", round(df.prf$var_explained/sum(df.prf$var_explained),3))
# View(df.prf)
for (i in 1:3){
  p <- ggplot(df.prf, aes(x="", y=var_explained, fill=traitPlus)) +
    geom_bar(stat="identity", width=1, color="black") +
    coord_polar("y", start=0) +
    theme_void() +
    # theme(legend.position="none") + 
    # geom_text(aes(y = ypos, label = var_explained ), color = "white", size=6) +
    scale_fill_brewer(palette = paste0("Set",i))
  print(p)
}
dev.off()

### FIG S8 ====
library(ggplot2)
load('dat_ulval/WGCNA_F_ulval_qreg_CP_reg2.RData')
all_f <- all; nf0 <- dim(all_f)[1]
load('dat_ulval/WGCNA_M_ulval_qreg_CP_reg2.RData')
all_m <- all; nm0 <- dim(all_m)[1]

modcorr <- merge(all_f[c(1,3)], all_m[c(1,3)], by="Geneid")
modcorr$count <- 1
modf <- unique(modcorr$module.x)
nf = length(modf)
modm <- unique(modcorr$module.y)
nm = length(modm)
datcorr <- as.data.frame(matrix(NA, ncol=3, nrow=nf*nm))
for (i in 1:nf){
  for (j in 1:nm){
    datcorr[(i-1)*nm+j,1] <- modf[i]
    datcorr[(i-1)*nm+j,2] <- modm[j]
    n1 <- sum(modcorr$module.x==modf[i])
    n2 <- sum(modcorr$module.y==modm[j])
    ncom <- sum(modcorr$module.x==modf[i] & modcorr$module.y==modm[j])
    datcorr[(i-1)*nm+j,3] <- ncom*min(nf0,nm0)/(n1*n2)
  } 
}
names(datcorr) <- c("female", "male", "odds")

pdf('dat_ulval/results/module_homology_FM.pdf')
ggplot(datcorr, aes(female, male)) + 
  geom_tile(aes(fill = odds), colour = "white") + 
  scale_fill_gradient(low = "white", high = "red", name="Fold Enrichment") +
  theme_grey(base_size = 8) +  
  labs(title='Correspondence between female and male PTSD modules', 
       x = "Female PTSD modules", y = "Male PTSD modules") +
  theme(axis.text.x = element_text(angle = 60, hjust = 1))
dev.off()

write.csv(datcorr[order(datcorr$odds, decreasing = T)[1:42],], 
          file="dat_ulval/results/module_homology_top.csv", row.names=F)

### FIG S9 ====
library(dplyr)
library(calibrate)
sa <- paste0(rep(c("F","M","A"),each=4), rep(c(11,25,24,9),3))
files <- paste0("dat_ulval/deseq_",sa,"_m3_ulval_keep3.RData")

pdf('results/volcano&heatmaps/volcano_MDD_deseq2.pdf')
for (k in 1:12){
  file = files[k]
  load(file)
  label <- strsplit(file,"_") %>% sapply("[",3) %>% gsub("r","",.)

  res <- sig01
  xmin <- -3.3  #min(res$log2FoldChange)-0.5
  xmax <- 3.3   #max(res$log2FoldChange)+0.5
  ymax <- 9.1 #max(-log10(res$pvalue))+0.5
  with(res,
       plot(log2FoldChange, -log10(pvalue), pch=20, main=label, col='grey', lwd=2,
            xlim=c(xmin,xmax), ylim=c(0,ymax), xlab="log2(FoldChange)"))
  with(subset(res, padj<.05 & log2FoldChange > 0),
       points(log2FoldChange, -log10(pvalue), pch=20, lwd=6, col="firebrick1"))#col="#FF6666"))
  with(subset(res, padj<.05 & log2FoldChange < 0),
       points(log2FoldChange, -log10(pvalue), pch=20, lwd=6, col="royalblue1"))#col="#006699"))
  if(min(res$padj)<.05){
    cric1 = max(res$pvalue[res$padj<=.05])
    cric2 = min(res$pvalue[res$padj>.05])
    abline(h=(-log10(cric1)-log10(cric2))/2, col="green", lty=2)
  }
}
dev.off()

### FIG S10 ====
rm(list=ls())
setwd('~/Documents/Projects/PTSD')
library(ggplot2)
library(dplyr)

pdf('results/wgcna&csea/module_PMcor_q_replot2.pdf', width=12, height=6)
labs = c("A","F","M"); labels <- c("Female+Male","Female","Male")
wfiles <- paste0('dat_ulval/WGCNA_', labs, '_ulval_qreg_CMP_reg2.RData')
for (k in 1:3){
  lab = labs[k]; label = labels[k]
  load(wfiles[k])
  modDiv <- data.frame(matrix(NA, nrow=dim(moduleTraitCor)[1], ncol=2))
  names(modDiv) <- c("module","PMcor")
  modDiv$module <- names(table(all$module))
  for (i in 1:dim(modDiv)[1]){
    mod <- modDiv$module[i]
    all.mod <- all[all$module==mod,]
    modDiv$PMcor[i] <- cor(all.mod$traitSignificance.GS.Dx.MDD, all.mod$traitSignificance.GS.Dx.PTSD)
  }
  modDiv$ngene <- sapply(modDiv$module, function(x) sum(all$module==x))
  modDiv$q <- .5
  N = dim(all)[1]
  K = 2000
  for (j in 1:dim(modDiv)[1]){
    n <- modDiv$ngene[j]
    p.dist <- sapply(1:K, function(x){
      sam <- sample(N,n)
      cor(all$traitSignificance.GS.Dx.MDD[sam], 
          all$traitSignificance.GS.Dx.PTSD[sam])
    })
    modDiv$q[j] <- sum(p.dist<modDiv$PMcor[j])/K
  }
  modDiv <- modDiv[order(modDiv$q, decreasing = T),]
  modDiv$col <- "black"
  col1 <- factor(modDiv$module, levels = modDiv$module)
  modDiv$col1 = col1
  p1 <- ggplot(modDiv, aes(x=col1, y=PMcor, fill=col1, col = col1)) + 
    geom_bar(position="dodge", stat="identity") +
    labs(title=label) +
    scale_fill_manual(values = modDiv$module) + 
    scale_color_manual(values = rep('black', dim(modDiv)[1])) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "none") 
  
  col2 <- factor(modDiv$module, levels = modDiv$module)
  col2 <- factor(col2, levels=rev(levels(col2)))
  modDiv$col2 = col2
  upper <- min(modDiv$q[modDiv$q>0.975])
  lower <- max(modDiv$q[modDiv$q<0.025])
  p2 <- ggplot(modDiv, aes(x=col2, y=q, fill=col2, col = col2)) + 
    geom_bar(position="dodge", stat="identity") +
    labs(title=label) +
    scale_fill_manual(values = modDiv$module) + 
    scale_color_manual(values = rep('black', dim(modDiv)[1])) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 60, hjust = 1),
          legend.position = "none") +
    geom_hline(yintercept = upper-0.01, col="blue") +
    geom_hline(yintercept = lower+0.01, col="red")
  p3 <- ggplot(data=modDiv, aes(x=col1, y=q, fill=col1, col=col1)) +
    geom_point(stat='identity', size=4)  + 
    scale_color_manual(values = modDiv$module) +
    xlab(NULL)+ylab("q") +
    theme_bw()+
    geom_hline(yintercept = upper-0.01, col="blue", lty=2) +
    geom_hline(yintercept = lower+0.01, col="red", lty=2) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1))
  print(p3)
}
dev.off()


### FIG S11a ====
library(ggplot2)
library(ggalt)
library(ggpubr)
rm(list=ls())
sex <- c('Female+Male','Female','Male')
cells <- c("N","A","O","M","E")
uns <- "U"
labs = c("A","F","M"); labels <- c("Female+Male","Female","Male")
wfiles <- paste0('dat_ulval/WGCNA_', labs, '_ulval_qreg_CMP_reg2.RData')
pdf('results/wgcna&csea/module_trait_cor_replot2.pdf')#, width=10, height=10)
for (k in 1:3){
  lab = labs[k]; label = labels[k]
  load(wfiles[k])
  
  moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
  moduleTraitPvalue$Dx.MDD.padj <- p.adjust(moduleTraitPvalue$Dx.MDD)
  moduleTraitPvalue$Dx.PTSD.padj <- p.adjust(moduleTraitPvalue$Dx.PTSD)
  mod.ord <- which(moduleTraitPvalue$Dx.MDD.padj <0.01 | moduleTraitPvalue$Dx.PTSD.padj<.01)# %>% as.vector
  moduleTraitCor <- as.data.frame(moduleTraitCor)
  moduleTraitCor <- moduleTraitCor[mod.ord,]
  moduleTraitPvalue <- moduleTraitPvalue[mod.ord,]
  moduleTraitCor$module <- rownames(to_plot)[mod.ord]
  colnames(to_plot) <- cells
  moduleTraitCor$cellmodule = paste0(uns,"_",moduleTraitCor$module)
  cnames <- colnames(to_plot)
  for (i in 1:dim(moduleTraitCor)[1]){
    if (sum(to_plot[mod.ord[i],]<0.01)){
      cell = sapply(to_plot[mod.ord[i],]<.01, isTRUE)
      moduleTraitCor$cellmodule[i] <- paste(c(cnames[cell], moduleTraitCor$module[i]), collapse="_")
    }
  }
  moduleTraitCor$err.ptsd = moduleTraitCor$Dx.PTSD[1]/qnorm(1 - moduleTraitPvalue$Dx.PTSD[1]/2)
  moduleTraitCor$err.mdd = moduleTraitCor$Dx.MDD[1]/qnorm(1 - moduleTraitPvalue$Dx.MDD[1]/2)
  col1 <- levels(factor(moduleTraitCor$module))
  moduleTraitCor$PTSD.Pvalue <- as.vector(moduleTraitPvalue$Dx.PTSD) %>% sprintf("%.3g",.)
  moduleTraitCor$MDD.Pvalue <- as.vector(moduleTraitPvalue$Dx.MDD) %>% sprintf("%.3g",.)
  moduleTraitCor$cellmodule <- factor(moduleTraitCor$cellmodule, levels=moduleTraitCor$cellmodule[order(moduleTraitCor$Dx.PTSD)])
  
  ### dotplot
  thres.m <- moduleTraitCor$err.mdd[1]; thres.p <- moduleTraitCor$err.ptsd[1]
  p1 <- ggplot(moduleTraitCor, aes(x=Dx.MDD, y=Dx.PTSD))+
    geom_point(shape=21, size=8,  color="black", fill=moduleTraitCor$module) +
    labs(title=paste0("Module-trait correlation, ", sex[k])) +
    geom_hline(yintercept = c(-1,1)*thres.p, lty=2, color="grey") +  
    geom_vline(xintercept = c(-1,1)*thres.m, lty=2, color="grey") +
    theme_classic2()
  p2 <- ggplot(moduleTraitCor, aes(x=Dx.MDD, y=Dx.PTSD))+
    geom_point(shape=21, size=5,  color="black", fill=moduleTraitCor$module) +
    labs(title=paste0("Module-trait correlation, ", sex[k])) +
    geom_text(data = moduleTraitCor, size = 4, colour = "black", fontface = "bold",
              aes(label = as.character(module)), hjust=.4, vjust=-1) +
    geom_hline(yintercept = c(-1,1)*thres.p, lty=2, color="grey") +  
    geom_vline(xintercept = c(-1,1)*thres.m, lty=2, color="grey") +
    theme_classic2()
  print(p1)
  print(p2)
}
dev.off()

### FIG S11b-d ====
library(gplots)
library(WGCNA)
files <- c('dat_ulval/WGCNA_A_ulval_qreg_CMP_reg2.RData',
           'dat_ulval/WGCNA_F_ulval_qreg_CMP_reg2.RData',
           'dat_ulval/WGCNA_M_ulval_qreg_CMP_reg2.RData')
sex <- c('Female+Male','Female','Male')
cells <- c("N","A","O","M","E")
uns <- "U"

filename = 'dat_ulval/results/module_csea_top_cell.pdf'
pdf(filename)
for (k in 1:3){
  load(files[k])
  moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
  moduleTraitPvalue$Dx.MDD.padj <- p.adjust(moduleTraitPvalue$Dx.MDD)
  moduleTraitPvalue$Dx.PTSD.padj <- p.adjust(moduleTraitPvalue$Dx.PTSD)
  mod.ord <- which(moduleTraitPvalue$Dx.MDD.padj <0.01 | moduleTraitPvalue$Dx.PTSD.padj<.01)# %>% as.vector
  moduleTraitCor <- as.data.frame(moduleTraitCor)
  moduleTraitCor <- moduleTraitCor[mod.ord,]
  moduleTraitPvalue <- moduleTraitPvalue[mod.ord,]
  moduleTraitCor$module <- rownames(to_plot)[mod.ord]
  colnames(to_plot) <- cells
  moduleTraitCor$cellmodule = paste0(uns,"_",moduleTraitCor$module)
  cnames <- colnames(to_plot)
  for (i in 1:dim(moduleTraitCor)[1]){
    if (sum(to_plot[mod.ord[i],]<0.01)){
      cell = sapply(to_plot[mod.ord[i],]<.01, isTRUE)
      moduleTraitCor$cellmodule[i] <- paste(c(cnames[cell], moduleTraitCor$module[i]), collapse="_")
    }
  }
  
  to_plot = to_plot[mod.ord,]
  rownames(to_plot) <- moduleTraitCor$cellmodule
  to_plot = to_plot[order(rownames(to_plot)),]
  idc <- apply(to_plot,1,min)<.01
  to_plot = to_plot[idc,]
  
  heatmap.2(-log10(to_plot),
            col=blueWhiteRed(1000,1)[500:1000],
            scale="none",
            trace="none",
            cexRow = 0.8,
            cexCol = .8, 
            density.info = "none",
            colsep=0:6,
            rowsep=0:sum(idc),
            sepcolor="grey",
            sepwidth=c(0.02,0.02),
            srtCol=45,
            offsetRow=0,
            offsetCol=-0.5,
            dendrogram="none",
            Rowv="none",
            Colv="none",
            key=T,
            key.xlab="-log10(P)",
            cellnote=signif(to_plot,1),
            notecex=.8,
            notecol="black",
            margins = c(12,10),
            main=paste0(sex[k], " CSEA"))
}
dev.off()

### END








