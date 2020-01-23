#### PTSD data analysis pipeline #### 2020-01-23

### DEG for sex- and region-specific data ----
library(DESeq2)

## for dat version 13
load('~/Documents/Projects/PTSD/data/dat13_ulval_keep3.RData')
cntdat <- dat13uv_keep[,-(1:6)]
demos <- meta13uv_keep
suff <- suff_keep

demos$RIN2 <- demos$RIN1^2
demos$Date[is.na(demos$Date)] <- "dat4"
demos$batch <- as.factor(demos$Date)
demos$Area <- demos$Region
demos$AgeDeath <- as.numeric(demos$AgeDeath)
demos$PMI <- as.numeric(demos$PMI)
rownames(demos) <- NULL

indices <- list(which(demos$Area==11 & demos$Sex=="F"),
                which(demos$Area==25 & demos$Sex=="F"),
                which(demos$Area==24 & demos$Sex=="F"),
                which(demos$Area==9  & demos$Sex=="F"),
                which(demos$Area==11 & demos$Sex=="M"),
                which(demos$Area==25 & demos$Sex=="M"),
                which(demos$Area==24 & demos$Sex=="M"),
                which(demos$Area==9  & demos$Sex=="M")
)
files <- paste(rep(c("F","M"),each=4), c(11,25,24,9), sep="")
alpha <- .5 ##by default .5

for (i in 1:8){
  index <- indices[[i]] 
  genes <- (rowSums(cntdat[,index]) > dim(cntdat)[2]*alpha)
  seldat <- cntdat[genes,index]
  dds <- DESeqDataSetFromMatrix(countData = seldat,
                                colData = demos[index,],
                                design= ~ PrimaryDx + AgeDeath + RIN1) ##model 3
  dds <- DESeq(dds)
  for (j in 1:3){
    if (j == 1) {
      res <- results(dds, contrast=c("PrimaryDx", "MDD", "Control"))
    } else if(j ==2) {
      res <- results(dds, contrast=c("PrimaryDx", "PTSD", "Control"))
    } else res <- results(dds, contrast=c("PrimaryDx", "PTSD", "MDD"))
    
    dd <- res@listData
    dd <- as.data.frame(dd)
    dd <- cbind(dd, data.frame(Geneid=suff$Geneid[genes], Genename=suff$Genename[genes]))
    dsort <- dd[order(dd$padj),]
    dsort$padj <- p.adjust(dsort$pvalue, method = "fdr")
    dsort <- dsort[!is.na(dsort$padj),]
    if (j == 1) {
      sig01 <- dsort
    } else if(j == 2){
      sig02 <- dsort
    } else sig12 <- dsort
  }
  fname <- paste0("deseq_", files[i], "_m3_ulval_keep3.RData")
  save(sig01, sig02, sig12, file=fname)
  
}

### pre-WGCNA normalization ----
## get gene annotation
getinfo <- c("ensembl_gene_id", "external_gene_id", "chromosome_name", "start_position", #"band", 
             "end_position", "strand", "gene_biotype", "percentage_gc_content")
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl", 
                host="feb2014.archive.ensembl.org") ## Gencode v19,
datProbes <- getBM(attributes=getinfo, filters=c("ensembl_gene_id"), values=dat13_u$Geneid, mart=mart)
save(datProbes, file="supplement/datProbes.RData")

library(biomaRt); library(WGCNA); library(corrplot); library(ggplot2); 
library(cqn); library(sva); library(limma); library(statmod); library(dplyr)

rm(list=ls())
load('supplement/datProbes.RData')
load('dat13_ulval_keep3.RData')
datExpr <- dat13r[,-(1:6)]
datMeta <- meta13r
datSuff <- dat13r[,1:6]

mat <- match(datSuff$Geneid, datProbes$ensembl_gene_id)
rownames(datProbes) = datProbes$ensembl_gene_id
datProbes$length = datProbes$end_position - datProbes$start_position
to_keep = !is.na(datProbes$length)
datProbes = datProbes[to_keep,]
gene_id <- match(datProbes$ensembl_gene_id, datSuff$Geneid)
datExpr <- datExpr[gene_id,]
rownames(datExpr) <- datProbes$ensembl_gene_id
datSuff <- datSuff[gene_id,]

sel_sam <- which(datMeta$PrimaryDx != "MDD" & datMeta$Sex=="F")
datMeta <- datMeta[sel_sam,]
datExpr <- datExpr[,sel_sam]

# ## FPKM filtering (no need if not qnormed)
# boxplot(datExpr,col=as.numeric(datMeta$Group), range=0)
# hist(apply(datExpr, 1, function(b) sum(b>.5)), breaks=50)
pres = apply(datExpr>1, 1, sum) 
to_keep = (pres > 0.5 * ncol(datExpr) & datProbes$gene_biotype=="protein_coding" & 
             datProbes$chromosome_name %in% c(1:22,"X","Y")) 
table(to_keep)
datExpr = datExpr[to_keep,]
datProbes = datProbes[to_keep,]

## CQN normalization
cqn.dat <- cqn(counts=datExpr, lengths=as.numeric(datProbes$length), 
               x=as.numeric(datProbes$percentage_gc_content),
               lengthMethod=c("smooth"), sqn=FALSE) ## Run cqn with specified depths and with no quantile normalization
cqn.dat <- cqn.dat$y + cqn.dat$offset ## Get the log2(Normalized FPKM) values
## y: The dependent value used in the systematic effect fit. Equal to log2 tranformed reads per millions.
## offset: The estimated offset. 
## sqn: This argument indicates whether the residuals from the systematic fit are (subset) quantile normalized. The default should only be changed by expert users.
datExpr.preCQN = datExpr
datExpr = cqn.dat

## remove outliers
normadj <- (0.5+0.5*bicor(datExpr))^2 ## Calculate connectivity
netsummary <- fundamentalNetworkConcepts(normadj)
ku <- netsummary$Connectivity
z.ku <- (ku-mean(ku))/sqrt(var(ku))
plot(1:length(z.ku),z.ku,col=datMeta$PrimaryDx,pch=19)
legend("bottomleft",legend = levels(factor(datMeta$PrimaryDx)), col = 1:3, pch=19, cex=.7)
abline(h=-4, lty=2)

outliers = (z.ku < -2)
table(outliers)
datExpr = datExpr[,!outliers]
datMeta = datMeta[!outliers,]


## Batch Correction
batch <- sapply(1:dim(datMeta)[1], function(x) substr(as.character(datMeta$Date)[x],1,4))
batch <- factor(batch)

datMeta$PrimaryDx <- factor(datMeta$PrimaryDx)
mod = model.matrix(~PrimaryDx, data=datMeta)
datExpr.combat = ComBat(dat = as.matrix(datExpr), batch = batch, mod = mod)
datExpr.cqn.preCombat = datExpr
datExpr = as.data.frame(datExpr.combat)

save(datExpr, datProbes, datMeta, 
     file="dat_ulval_F_qnorm_CMP.RData")



### WGCNA for sex-specific data ----
rm(list=ls())
load('dat_ulval_F_qnorm_CMP.RData')
dat <- t(datExpr)
demos <- datMeta
# id.sel <- which(demos$Sex=="F")
# demos <- demos[id.sel,]
# dat <- dat[id.sel,]

## library
library(magrittr)
library(stats)
library(WGCNA)
library(flashClust)
library(biomaRt)
options(stringsAsFactors = FALSE, digits = 3)
library(corrplot)
library(ggplot2)
library(cqn)
library(sva)
library(limma)
library(statmod)
theme_update(plot.title = element_text(hjust = 0.5))
library(dplyr)
library(pSI)
library(pSI.data)
library(gplots)


##
demos$Area <- demos$Region
demos$batch <- demos$Date
datMeta <- demos[c("AgeDeath", "PMI", "Race", "Smoking", "batch", "RIN1")]
datMeta$Dx.MDD <- NA; datMeta$Dx.MDD[demos$PrimaryDx=="MDD"] <- 1; datMeta$Dx.MDD[demos$PrimaryDx=="Control"] <- 0
datMeta$Dx.PTSD <- NA; datMeta$Dx.PTSD[demos$PrimaryDx=="PTSD"] <- 1; datMeta$Dx.PTSD[demos$PrimaryDx=="Control"] <- 0
datMeta$PrimaryDx <- as.numeric(demos$PrimaryDx)
datMeta$Race <- as.character(datMeta$Race); datMeta$Race[datMeta$Race=="W"] <- 1; datMeta$Race[datMeta$Race=="B"] <- -1; datMeta$Race[datMeta$Race!=1&-1] <- 0 
datMeta$Race <- as.numeric(datMeta$Race)
datMeta$Smoking <- as.character(datMeta$Smoking); datMeta$Smoking[datMeta$Smoking=="Y"] <- 1; datMeta$Smoking[datMeta$Smoking=="N"] <- -1; datMeta$Smoking[datMeta$Smoking=="U"] <- 0
datMeta$Smoking <- as.numeric(datMeta$Smoking)
datMeta$BA11 <- as.numeric(demos$Area==11)
datMeta$BA24 <- as.numeric(demos$Area==24)
datMeta$BA25 <- as.numeric(demos$Area==25)
datMeta$BA9 <- as.numeric(demos$Area==9)
datMeta$batch <- as.numeric(datMeta$batch)
#datMeta$Sex <- as.numeric(demos$Sex)
traitColors <- numbers2colors(datMeta, signed = FALSE)


##** regression model 2: reg2**
lm.cov <- lm(dat ~ AgeDeath + RIN1 + PMI, data=datMeta) ##for sex specific, not including Race/Area
dat.ori <- dat
dat <- dat - as.matrix(datMeta[c('AgeDeath', 'RIN1', 'PMI')]) %*% lm.cov$coefficients[2:4,]


##
enableWGCNAThreads()
powers <- c(c(1:10), seq(from = 12, to=20, by=2))
sft <- pickSoftThreshold(dat, powerVector = powers, verbose = 3)
#sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h = 0.9, col = "red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
print(sft$fitIndices[,5])


## module construction
R2 <- -sign(sft$fitIndices[,3])*sft$fitIndices[,2]
power <- min(which(R2>.85))
#power <- 5
net <- blockwiseModules(dat, power = power,
                        TOMType = "signed", minModuleSize = 20,
                        reassignThreshold = 0, mergeCutHeight = 0.1,
                        numericLabels = TRUE, pamRespectsDendro = F,
                        saveTOMs = F, 
                        saveTOMFileBase = "catTOM",
                        verbose = 3)

#plotDendroAndColors(net$dendrograms[[1]], net$colors[net$blockGenes[[1]]], 
#    main = "Single block gene dendrogram and module colors",
#    dendroLabels = FALSE, hang = 0.03,
#    addGuide = TRUE, guideHang = 0.05)
#plotDendroAndColors(net$dendrograms[[2]], net$colors[net$blockGenes[[2]]], 
#    main = "Single block gene dendrogram and module colors",
#    dendroLabels = FALSE, hang = 0.03,
#    addGuide = TRUE, guideHang = 0.05)
#plotDendroAndColors(net$dendrograms[[3]], net$colors[net$blockGenes[[3]]], 
#    main = "Single block gene dendrogram and module colors",
#    dendroLabels = FALSE, hang = 0.03,
#    addGuide = TRUE, guideHang = 0.05)


## module-trait cor
moduleColors <- labels2colors(net$colors)
nGenes = ncol(dat);
MEs0 = moduleEigengenes(dat, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
nSamples <- dim(datMeta)[1]
moduleTraitCor = cor(MEs, datMeta, use = "p");# Pearson correlation
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
textMatrix = paste(signif(moduleTraitCor, 2), " (",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(5, 10, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datMeta),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(dat, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
# geneTraitSignificance = as.data.frame(cor(dat, mydataP, use = "p"));
geneTraitSignificance = as.data.frame(cor(dat, datMeta, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(datMeta), sep="");
names(GSPvalue) = paste("p.GS.", names(datMeta), sep="");
all <- data.frame(Geneid=datProbes$ensembl_gene_id,
                  Genename=datProbes$external_gene_id,
                  module=moduleColors)
all <- cbind(all, traitSignificance=geneTraitSignificance)
mm <- sapply(1:dim(all)[1], function(x){geneModuleMembership[x, match(moduleColors[x], modNames)]})
all$DxmoduleMembership <- mm


## CSEA
load('zhang.pSIout.RData')
colors = as.vector(all$module)
cell.p.zhang = matrix(NA, length(unique(colors)), 5);  rownames(cell.p.zhang) = unique(colors)
colnames(cell.p.zhang) = colnames(pSI.output)
for(mod in unique(colors)) {
  f = fisher.iteration(pSI.output, all$Geneid[all$module==mod], p.adjust = F)
  cell.p.zhang[mod,] = f$`0.05 - nominal`
}

# modSig <- rownames(moduleTraitPvalue)[moduleTraitPvalue[,1]<.05]
# modSig <- gsub('ME','',modSig)
modSig <- modNames
names(MEs) <- gsub('ME', '', names(MEs))
cell.p.zhang.fdr = p.adjust(cell.p.zhang,"fdr")
dim(cell.p.zhang.fdr) = dim(cell.p.zhang); dimnames(cell.p.zhang.fdr) = dimnames(cell.p.zhang);
indSig <- match(modSig, rownames(cell.p.zhang.fdr))
to_plot = cell.p.zhang.fdr[indSig,]
dendro.col = dendro.col.zhang
denro.row= as.dendrogram(hclust(as.dist(1-bicor(MEs[modSig])), method="average"))

heatmap.2(-log10(to_plot+10^(-100)),
          col=blueWhiteRed(1000,1)[500:1000],
          scale="none",
          trace="none",
          cexRow = 0.8,
          cexCol = .8, 
          density.info = "none",
          colsep=0:7,
          rowsep=0:length(modSig),
          sepcolor="grey",
          sepwidth=c(0.02,0.02),
          srtCol=45,
          offsetRow=0,
          offsetCol=-0.5,
          Rowv=denro.row,
          Colv=dendro.col,
          key=T,
          key.xlab="-log10(P)",
          cellnote=signif(to_plot,1),
          notecex=.8,
          notecol="black",
          #           main=paste0(sdlab," Enrichment")) ###for S x Dx
          #           main=paste0(salab, " Enrichment")) ###for S x Area
          main='Enrichment')


##save
filesave = "WGCNA_F_ulval_qreg_CMP_reg2.RData"
save(dat, datProbes, datMeta, all, to_plot, denro.row, dendro.col, moduleTraitCor, moduleTraitPvalue, MEs, textMatrix, 
     file=filesave)