## Key driver analysis ##
#####################################
# Load packages and useful functions
library('grDevices')
library('WGCNA')
library('bnlearn')
library('KDA')
library('qgraph')
library('plyr')
library('pracma')
library('ggplot2')
library('YaleToolkit')
library('igraph')
library('circlize')
library('tidyr')
library('dplyr')
library('grDevices')
library('ggplotify')
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# Load files and prepare paths for output, change it if needed
args <- c("/home/jw2372/work/PTSD/dat_ulval/dat13_ulval_keep3.RData", "/home/jw2372/work/PTSD/dat_ulval/", "deseq_[AFM][0-9]+_m3_ulval_keep3.RData", 
          "/home/jw2372/work/PTSD/dat_ulval/", "WGCNA_[AFM]+_ulval_qreg_CP_reg2.RData", "/SAY/standard2/IAC-CC0937-BiostatYSPH/dj333/.Project/PTSD/PTSD/Results/Figures/", "0826")
#####################################
# Change it if you want to save files in your folder
rootdir <- '/SAY/standard2/IAC-CC0937-BiostatYSPH/dj333/.Project/PTSD/PTSD/'
data_file <- args[1]
deg_file_path <- args[2]
deg_file_pattern <- args[3]
deseq_files <- dir(path = deg_file_path, pattern = deg_file_pattern)
area_by_sex <- gsub(".+_([AFM][0-9]+).+", "\\1", deseq_files)
tmp_sex_codes <- gsub("[^AFM]", "", area_by_sex)
sex_codes <- unique(tmp_sex_codes)
area_codes <- sort(unique(as.numeric(gsub("[AFM]", "", area_by_sex))))

wgcna_file_path <- args[4]
wgcna_file_pattern <- args[5]
wgcna_files <- dir(path = wgcna_file_path, pattern = wgcna_file_pattern)
wgcna_sex_codes <- gsub(".+_([AFM]).+", "\\1", wgcna_files)
wgcna_files <- paste0(wgcna_file_path, "/", wgcna_files)

sex_codes <- intersect(sex_codes, wgcna_sex_codes)

area_by_sex <- area_by_sex[tmp_sex_codes %in% sex_codes]
deseq_files <- deseq_files[tmp_sex_codes %in% sex_codes]
if(length(area_by_sex) != length(sex_codes) * length(area_codes)){
    warning("Each sex code should have the same list of areas!")
}

deseq_files <- paste0(deg_file_path, "/", deseq_files)
output_file_path <- args[6]
outfile_suffix <- args[7]
# Create paths for output
mainDir <- paste0(output_file_path, "/", outfile_suffix)
dir.create(mainDir, recursive = T, showWarnings = F)
Figures_path <- file.path(mainDir, "Figures")
Tables_path <- file.path(mainDir, "Tables")
RData_path <- file.path(mainDir, "RData")
dir.create(Figures_path, showWarnings = F)
dir.create(Tables_path, showWarnings = F)
dir.create(RData_path, showWarnings = F)
#####################################
##### If your files do not match the requirements, please make the adjustments here #####
fpkm <- fpkm13uv_keep
#fpkm <- cbind(suff_keep[, 1:6], fpkm)
demos <- meta13uv_keep
#demos$Sample <- demos$Samples
#####################################
fpkmdat <- as.data.frame(fpkm[, as.character(demos$Sample)])
#fpkmdat <- t(apply(fpkmdat, 1, function(x){return((x - mean(x)) / sd(x))}))
fpkmdat <- cbind(fpkm[, 1:6], fpkmdat)
rownames(demos) <- demos$Sample
demos$PrimaryDx <- as.character(demos$PrimaryDx)

# Set colors: light colors for down regulation and dark colors for up regulation
# Mad
colors <- matrix(c('pink', 'white', 'red', # Color for 4
                   'lightblue', 'white', 'blue', # Color for 9
                   'yellow', 'white', 'orange', # Color for 11
                   'lightgreen', 'white', 'green', # Color for 24
                   'azure3', 'white', 'azure4' ), ncol= 3, byrow = T) # Color for 25

# Set sex code
sexcode <- list()
sexcode$F <- "females"
sexcode$M <- "males"
sexcode$A <- 'all'

sexcode1 <- list()
sexcode1$F <- "Female"
sexcode1$M <- "Male"
sexcode1$A <- 'All'

# Set potential groups
g <- c(paste0('MDD.FM', sort(area_codes)), paste0('PTSD.FM', sort(area_codes)))

# Set vertices size
v.sizes <- linspace(x1 = 4, x2 = 20, 40)
#####################################
# Aggregate DEG information by sex group.
all0102 <- data.frame(Geneid = fpkmdat$Geneid)
for(i in 1:length(area_by_sex)){
    code <- area_by_sex[i]
    load(deseq_files[i])
    sig01[, paste0('MDD.', code)] <- as.numeric(sig01$padj < 0.1) * sign(sig01$log2FoldChange)
    sig02[, paste0('PTSD.', code)] <- as.numeric(sig02$padj < 0.1) * sign(sig02$log2FoldChange)
    sig0102 <- join(x = sig01, y = sig02, by = 'Geneid', 'left', 'first')
    all0102 <- join(all0102, sig0102[, c('Geneid', paste0('MDD.', code), paste0('PTSD.', code))], by = 'Geneid', 'left', 'first')
}

case_by_sex <- unique(gsub("[0-9]+$", "", colnames(all0102)[-1]))
for(case in case_by_sex){
       all0102[, case] <- apply(all0102[, grep(case, colnames(all0102))] != 0, 1, any)
}

# Area names
area_name_list <- rep("unknown", 25)
area_name_list[9] <- 'dlPFC'
area_name_list[11] <- 'OFC'
area_name_list[24] <- 'dACC'
area_name_list[25] <- 'sgPFC'
#####################################
# Global variable "datMeta0" for demographic data
datMeta0 <- demos
# Main function for analysis
group <- 'A'
modules_mult_scale <- 4 
sizes <- c(1, 1, 1, 1, 1, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8)
KDA_v1 <- function(group){
	## Use aggregated DEG information "all0102", 
	## demographic information "datMeta0", and expression
	## profile "fpkmdat"
	## The function to output:
	1. KDA Results
	2. Module DEG enrichment results
	3. Organized R tables for further visualization (circos)
	4. PDF files for network visualization with DEG and KDA
    groups <- gsub("FM", group, g)
    code <- group
    load(wgcna_files[grep(gsub("AFM", code, "_AFM"), wgcna_files)])
    # Collect samples with their information  
    demos <- datMeta0
    if(code == 'A'){
        demos$Sex <- group
    } else{
        samples <- as.character(demos$Sample[demos$Sex == group])
        demos <- datMeta0[samples, ]
    }
    module_eigen <- matrix(0, nrow=0, ncol = 7)
    ########################################
    # Get information of module significance
    rownames(moduleTraitPvalue) <- gsub('ME', '', rownames(moduleTraitPvalue))
    moduleTraitPvalue <- as.data.frame(moduleTraitPvalue)
    moduleTraitPvalue$Dx.MDD.fdr <- p.adjust(moduleTraitPvalue$Dx.MDD, 'fdr')
    moduleTraitPvalue$Dx.PTSD.fdr <- p.adjust(moduleTraitPvalue$Dx.PTSD, 'fdr')
    modules_csea_fdr <- to_plot
    ########################################
    dend <- denro.row
    minfo <- as.data.frame(all[, c('Geneid', 'module', 'Genename')], stringsAsFactors = F)
    colnames(minfo) <- c('Geneid', 'module', 'Genename')
    minfo <- join(minfo, all0102[, c(1, grep(paste0('\\.', group), colnames(all0102)))], by = 'Geneid', 'left', 'first')
    # Print the number of DEGs 
    print(sprintf('For this dataset, there are %d genes differentially expressed for MDD and %d genes for PTSD', 
           sum(minfo[, paste0('MDD.', group)], na.rm = T), sum(minfo[, paste0('PTSD.', group)], na.rm = T)))
    minfo$Genename <- make.names(minfo$Genename, unique=T)
    # Get module names
    modules <- names(table(minfo$module))
    moduleTraitPvalue <- moduleTraitPvalue[modules, ]
    modules_csea_fdr <- modules_csea_fdr[modules, ]
    to_plot <- to_plot[modules, ]
    modules.enrichment <- matrix(NA, ncol = ncol(minfo)-3, nrow = length(modules))
    colnames(modules.enrichment) <- colnames(minfo)[-c(1:3)]
    modules.enrichment.up <- modules.enrichment
    modules.enrichment.down <- modules.enrichment
    
    print("here")
        
    # Run fisher exact test on each module for DEGs (up/down/up or down) at each region
    # calculate MDC for each module
    MDC <- rep(1, length(modules))
    MDC1 <- rep(1, length(modules))
    cons <- matrix(NA, nrow = length(modules), ncol = 3)
    for(k in 1:length(modules)){
        print(k)
        m <- modules[k]
        module <- as.numeric(minfo$module == m)
        # Extract gene expression info from samples in Control
        datExpr <- join(minfo[minfo$module == m, 'Geneid', drop = F], fpkmdat[, c('Geneid', demos$Sample[which(demos$PrimaryDx == 'Control' & demos$Sex == group)])], 
                        by = 'Geneid', 'left', 'first')
        datET0 <- datExpr[, -grep('Geneid', colnames(datExpr))]
        plot_genes <- rownames(datET0)
        
        # Extract gene expression info from samples in PTSD
        datExpr <- join(minfo[minfo$module == m, 'Geneid', drop = F], fpkmdat[, c('Geneid', demos$Sample[which(demos$PrimaryDx == 'PTSD' & demos$Sex == group)])], 
                        by = 'Geneid', 'left', 'first')
        datET1 <- datExpr[, -grep('Geneid', colnames(datExpr))]
    
        # Extract gene expression info from samples in MDD
        datExpr <- join(minfo[minfo$module == m, 'Geneid', drop = F], fpkmdat[, c('Geneid', demos$Sample[which(demos$PrimaryDx == 'MDD' & demos$Sex == group)])], 
                        by = 'Geneid', 'left', 'first')
        datET2 <- datExpr[, -grep('Geneid', colnames(datExpr))]
    
        
        
            
        for(item in colnames(minfo)[-c(1:3)]){
            MDDorPTSD <- factor(as.numeric(minfo[, item] != 0), levels = c(0, 1))
            modules.enrichment[k, item] <- fisher.test(table(module, MDDorPTSD), alternative = 'greater')$p.value
            if(length(grep('^[A-Z]+\\.[A-Z]+$', item)) == 1){
                MDDorPTSD.up <- factor(as.numeric(apply(minfo[, grep(paste0(item, '[0-9]+'), colnames(minfo))] == 1, 1, any)), levels = c(0, 1))
                MDDorPTSD.down <- factor(as.numeric(apply(minfo[, grep(paste0(item, '[0-9]+'), colnames(minfo))] == -1, 1, any)), levels = c(0, 1))
            }
            else{
                MDDorPTSD.up <- factor(as.numeric(minfo[, item] == 1), levels = c(0, 1))
                MDDorPTSD.down <- factor(as.numeric(minfo[, item] == -1), levels = c(0, 1))
            }
            modules.enrichment.up[k, item] <- fisher.test(table(module, MDDorPTSD.up), alternative = 'greater')$p.value
            modules.enrichment.down[k, item] <- fisher.test(table(module, MDDorPTSD.down), alternative = 'greater')$p.value
        }
    }
    colnames(modules.enrichment.up) <- paste0(colnames(modules.enrichment), '.up')
    colnames(modules.enrichment.down) <- paste0(colnames(modules.enrichment), '.down')
    modules.enrichment <- cbind(modules.enrichment, modules.enrichment.up, modules.enrichment.down)
    rownames(modules.enrichment) <- modules
    modules.enrichment.p <- modules.enrichment
    colnames(modules.enrichment.p) <- paste0(colnames(modules.enrichment), ".p")
    modules.enrichment.p <- as.data.frame(modules.enrichment.p)

    for(i in 1:ncol(modules.enrichment.p)){
        modules.enrichment[, i] <- p.adjust(modules.enrichment.p[, i], method = 'fdr')
    }
    colnames(modules.enrichment) <- paste0(colnames(modules.enrichment), ".fdr")     
    print(paste0(sprintf('For this dataset, there are %d modules with PTSD DE genes enriched: ', sum(modules.enrichment[, paste0('PTSD.', group, '.fdr', collapse = '')] < 0.05)),
           paste0(modules[which(modules.enrichment[, paste0('PTSD.', group, '.fdr', collapse = '')] < 0.05)], collapse = ' , '), '.'))
    
    # KDA analysis for modules with significant enrichment
    key.drivers <- list()
    pdf(file = paste0(Figures_path, '/Modules_', group, '_', outfile_suffix, '_KDA.pdf'), width = 10, height = 10)
    csea_sig <- apply(modules_csea_fdr < 0.05, 1, any)
    trait_sig <- (moduleTraitPvalue$Dx.MDD.fdr < 0.05) | (moduleTraitPvalue$Dx.MDD.fdr < 0.05)
    PTSD_enriched_sig <- modules.enrichment[, paste0(c('PTSD.', group, '.fdr'), collapse = '')] < 0.05
    MDD_enriched_sig <- modules.enrichment[, paste0(c('PTSD.', group, '.fdr'), collapse = '')] < 0.05
    any_sig <- csea_sig | trait_sig | PTSD_enriched_sig | MDD_enriched_sig
    if(sum(modules.enrichment[, paste0(c('PTSD.', group, '.fdr'), collapse = '')] < 0.05) >= 1){
    # Significant ones only
    #    for(m in modules[modules.enrichment[, paste0(c('PTSD.', group, '.fdr'), collapse = '')] < 0.05]){
    # All modules ranked by orders
        for(m in modules[order(modules.enrichment[, paste0(c('PTSD.', group, '.fdr'), collapse = '')])]){         
            tryCatch({     
            # Gather expression profile and DEG information for genes in this module
            DEm <- minfo[which(minfo$module == m & minfo[, paste0('PTSD.', group)]), 'Genename',drop=F]
            datExpr <- join(minfo[minfo$module == m & !is.na(minfo[, paste0('PTSD.', group)]), c('Geneid', 'Genename'), drop = F], fpkmdat[, c('Geneid', demos$Sample[which(demos$Sex == group)])], 
                            by = 'Geneid', 'left', 'first')
            datET <- collapseRows(datET = datExpr[, -grep('Gene', colnames(datExpr))], rowGroup = datExpr[, 'Genename'], rowID = rownames(datExpr), connectivityBasedCollapsing = T)$datETcollapsed
            plot_genes <- rownames(datET)
            print(paste0('Module ', m, ' has ', dim(datET)[1], ' genes and ', 
                         table(as.numeric(minfo$module == m), minfo[, paste0('PTSD.', group)])[2, 2], 
                         ' PTSD DE genes (FDR ', sprintf("%.2e)", modules.enrichment[m, paste0('PTSD.', group, '.fdr', collapse = '')])))
            
            # Infer the graph using aracne
            #aracne.results <- aracne(as.data.frame(t(datET)))
            datET.pca <- prcomp(t(datET), center = T,scale. = T)
            sample_tmp <- demos$Sample[which(demos$Sex == group)]
            pred.MDD <- predict(datET.pca, newdata=t(datET[, demos[sample_tmp, 'PrimaryDx'] == "MDD"]))
            pred.Control <- predict(datET.pca, newdata=t(datET[, demos[sample_tmp, 'PrimaryDx'] == "Control"]))
            pred.PTSD <- predict(datET.pca, newdata=t(datET[, demos[sample_tmp, 'PrimaryDx'] == "PTSD"]))
            module_eigen <- rbind(module_eigen, c(m, mean(pred.Control[, 1]), sd(pred.Control[, 1]), mean(pred.MDD[, 1]), sd(pred.MDD[, 1]), mean(pred.PTSD[, 1]), sd(pred.PTSD[, 1])))
            aracne.results <- aracne(as.data.frame(t(datET)))
            connection_m <- aracne.results$arcs
            connection_m <- connection_m[connection_m[, 2] < connection_m[, 1], ,drop=F]
            # Run key driver analysis based on the graph
            if(dim(datExpr)[1] >= 10){
                KDA_results <- KDA::keyDriverAnalysis(connection_m, signature = as.character(DEm[,1]), directed = F, nlayer_search=3, expanded_network_as_signature=T)
                if(is.null(KDA_results)){
                    cols <- c('keydrivers','is_signature','hits','downstream','signature_in_subnetwork','subnetwork_size','signature_in_network','network_size','signature','optimal_layer','fold_change_whole','pvalue_whole','fold_change_subnet','pvalue_subnet','pvalue_corrected_subnet','keydriver','Module')
                    key.driver <- as.data.frame(matrix(NA, nrow = 0, ncol = length(cols)))
                    colnames(key.driver) <- cols
                }
                else{
                    key.driver <- as.data.frame(KDA_results$keydrivers, stringsAsFactors = F)
                    key.driver$Module <- m
                    key.driver <- key.driver[key.driver$is_signature == 1, ]
                    key.drivers[[length(key.drivers)+1]] <- key.driver
                }
            }
            if(dim(datET)[1] > 50){
                DEonly <- apply(aracne.results$arcs, 1, function(x) (x[1] %in% DEm[, 1]) | (x[2] %in% DEm[, 1]))
                connection_m <- aracne.results$arcs[DEonly,,drop=F]
                connection_m <- connection_m[connection_m[, 2] < connection_m[, 1], ,drop=F]
                plot_genes <- unique(c(aracne.results$arcs[DEonly, 'from'], aracne.results$arcs[DEonly, 'to']))
            }
            if(length(plot_genes) > 100){
                DEonly <- apply(connection_m, 1, function(x) (x[1] %in% c(key.driver$keydrivers, DEm[, 1])) & (x[2] %in% c(key.driver$keydrivers, DEm[, 1])))
                connection_m <- connection_m[DEonly, ,drop=F]
                plot_genes <- unique(as.character(connection_m))
            }
            # Get the annotation
            csea_anno <- paste0('CSEA: ', paste0(colnames(modules_csea_fdr)[modules_csea_fdr[m, ] < 0.05], collapse = ', '))
            trait_anno <- paste0('Trait: ', paste0(c('MDD', 'PTSD')[moduleTraitPvalue[m, c('Dx.MDD.fdr', 'Dx.PTSD.fdr')] < 0.05], collapse = ', '))                     
            # Plot the network
            graph_m <- graph_from_edgelist(connection_m, directed = F)
            vertex.frame.color <- 'black'
            if(dim(datExpr)[1] >= 10 & dim(connection_m)[1] > 0){
                    tryCatch({
                vertex.frame.color <- lapply(V(graph_m)$name, function(x) {
                                        if(x %in% key.driver$keydrivers){
                                            return('magenta')
                                        }
                                        else return('black')
                                    }
                                    )
                vertex.frame.color <- do.call(c, vertex.frame.color)
                vertex.size.DE <- lapply(V(graph_m)$name, function(x) {
                                        if(x %in% DEm[,1]){
                                            return(2 * modules_mult_scale)
                                        }
                                        else return(0.5)                                    
                               })
                vertex.size.DE <- do.call(c, vertex.size.DE)
                vertex.label.dist.DE <- lapply(V(graph_m)$name, function(x) {
                                        if(x %in% DEm[,1]){
                                            return((0.3) * modules_mult_scale)
                                        }
                                        else return(0.1)                                    
                               })
                vertex.label.dist.DE <- do.call(c, vertex.label.dist.DE)                        
                vertex.label.cex.DE <- lapply(V(graph_m)$name, function(x) {
                                        if(x %in% DEm[,1]){
                                            return((0.2) * modules_mult_scale * sizes[nchar(x)])
                                        }
                                        else return((0.2) * modules_mult_scale * sizes[nchar(x)])                                    
                               })
                vertex.label.cex.DE <- do.call(c, vertex.label.cex.DE)                          
                
                vertex.label.dist=(0.3) * modules_mult_scale
                V(graph_m)$label.cex = vertex.label.cex.DE
                #V(graph_m)$vertex.size = vertex.size.DE
                pie.values <- lapply(V(graph_m)$name, function(x) rep(1, length(names(which(table(demos$Area) > 0)))))
                pie.color <- lapply(V(graph_m)$name, function(x) {
                                tmp <- as.numeric(minfo[minfo$Genename == x, grep('PTSD.[A-Z]+[0-9]+', colnames(minfo))] + 2)
                                return(colors[as.matrix(cbind(1:length(tmp), tmp))])
                                })
                pvalue <- sprintf("%.2e", modules.enrichment[m, paste0('PTSD.', group, '.fdr', collapse = '')])
                e <- get.edgelist(graph_m, names = F)
                l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(graph_m), area=8*(vcount(graph_m)^2),repulse.rad=(vcount(graph_m)^3.1))
                #l <- layout_with_fr(graph_m)
                plot(graph_m, layout = l, vertex.label.dist=vertex.label.dist.DE, vertex.label.degree=pi/2, vertex.shape='pie', vertex.frame.color = vertex.frame.color,
                     vertex.pie = pie.values, vertex.pie.color=pie.color, vertex.size = vertex.size.DE, vertex.label.font=2) +
                     #vertex.pie = pie.values, vertex.pie.color=pie.color, vertex.size = v.sizes[igraph::degree(graph_m)]) +
                                     title(main = paste0('Module ', m, ' of all ', sexcode[[group]], " (p-value=", pvalue, ") \n", csea_anno, "\n", trait_anno)) }, warning = function(w) {print(paste("Plot failure:", m))},  error = function(e) {print(e)})
            }
        }, warning = function(w) {print(paste("Data failure:", m))},  error = function(e) {print(e)})                                
        }
    }
    dev.off()
    
    # Prepare data for circos plot
    # Select module enrichment info with respect to PTSD for all regions                              
    tmp <- paste0(as.character(t(cbind(colnames(modules.enrichment.up), colnames(modules.enrichment.down)))), '.fdr')                     
    colfunc<-colorRampPalette(c("yellow","red"))
    tmp <- tmp[grep('PTSD', tmp)]
    
    ##### Results of CSEA                             
    m2plot <- to_plot
    colfunc<-colorRampPalette(c("yellow","red"))
    m2plot <- -log10(m2plot)
    cuts <- c(1.3, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    m2plot1 <- m2plot
    i <- 1
    m2plot1[m2plot < cuts[i]] <- 'azure3'
    for(i in 1:9){
        m2plot1[m2plot >= cuts[i] & m2plot < cuts[i+1]] <- colfunc(10)[i]
    }
    i <- 10
    m2plot1[m2plot >= cuts[i]] <- colfunc(10)[i]
    m2plot_csea <- m2plot1  
    
    ######
    
    m2plot <- modules.enrichment[, tmp]
    m2plot <- -log10(m2plot)
    cuts <- c(1.3, 2, 3, 4, 5, 6, 7, 8, 9, 10)
    m2plot1 <- m2plot
    i <- 1
    m2plot1[m2plot < cuts[i]] <- 'azure3'
    for(i in 1:9){
        m2plot1[m2plot >= cuts[i] & m2plot < cuts[i+1]] <- colfunc(10)[i]
    }
    i <- 10
    m2plot1[m2plot >= cuts[i]] <- colfunc(10)[i]
    m2plot1 <- as.data.frame(m2plot1, stringsAsFactors = F)
    m2plot2 <- m2plot1
    m2plot1$MDC <- (MDC - min(MDC)) / (max(MDC) - min(MDC))
    m2plot1$MDC1 <- (MDC1 - min(MDC1)) / (max(MDC1) - min(MDC1))
    m2plot1$MDC2 <- (MDC2 - min(MDC2)) / (max(MDC2) - min(MDC2))
    
    # Prepare  legends
    m2plot2[1:11, 3] <- c('azure3', colfunc(10))
    m2plot2$module <- rownames(m2plot1)
    m2plot2 <- gather(m2plot2, key, value, -module)
    m2plot2$value <- factor(m2plot2$value, levels = rev(c('azure3', colfunc(10))))
    plot2 <- ggplot(data = m2plot2, aes(x = module, y = key)) + geom_tile(aes(fill = value)) + 
            scale_fill_manual('', values = rev(c('azure3', colfunc(10))), 
                              labels = rev(c('NS', '1.3', '2', '3', '4', '5', '6', '7', '8', '9', '>=10'))) + 
    theme(legend.text=element_text(size=10))
    scale_bar <- g_legend(plot2)
     
    cell_names<- c('Neuron', 'Astrocyte', 'Oligodendrocyte', 'Microglia', 'Endothelial')
    A <- as.data.frame(matrix(sample(1:(((length(names(which(table(demos$Area) > 0)))) + 1)*2 + 5), replace = T, size = 100), nrow = 10))
    A$module <- paste0('V', 1:10)
    A <- gather(A, key, value, -module)
    A$value <- paste0('V', A$value)
    plot3 <- ggplot(data = A, aes(x = module, y = key)) + geom_tile(aes(fill = value), color = "white") + 
            scale_fill_manual('', values = c(rep('azure3', (((length(names(which(table(demos$Area) > 0)))) + 1)*2 + 5))), 
                              labels = rev(c(cell_names, as.character(t(cbind(paste0('', area_name_list[sort(as.numeric(names(which(table(demos$Area) > 0))))],' up'), 
                                                                  paste0('', area_name_list[sort(as.numeric(names(which(table(demos$Area) > 0))))],' down')))),
                                                                  'ALL up', 'ALL down'))) + theme(legend.text=element_text(size=10))
    legend_bar <- g_legend(plot3)  
    write.table(x = do.call(rbind, key.drivers), file = paste0(Tables_path, '/KDA_', group, '_', outfile_suffix, '_KDA.csv'), quote=F, sep = ",", row.names = F)                        
    return(list(enrichment = modules.enrichment, csea = modules_csea_fdr,
                plots = list(plot_mat = m2plot1, plot_mat_csea = m2plot_csea, dend = dend, cons = cons, scale_bar = scale_bar, legend_bar = legend_bar),
                key.drivers = key.drivers, module_eigen = module_eigen))                            
}
#####################################




## Compare with Gandal's Meta Analysis ##
### The code below is directly adapted from github repo "Shared-molecular-neuropathology-across-major-psychiatric-disorders-parallels-polygenic-overlap", (Gandal, M. J., Haney, J. R., Parikshak, N. N., Leppa, V., Ramaswami, G., Hartl, C., ... & Liu, C. (2018). Shared molecular neuropathology across major psychiatric disorders parallels polygenic overlap. Science, 359(6376), 693-697.) 
# Compare our DEG with gandall's (MDD/PTSD)
PLOT_PDF=F
areas <- area_by_sex
plots <- list()
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
    all_genes = intersect(all_genes, MDD_meta$Geneid)
    allmeta = matrix(NA,nrow=length(all_genes), 8)
    allmeta[,1] = asd_meta$beta[match(all_genes, rownames(asd_meta))]
    allmeta[,2] = scz_meta$beta[match(all_genes, rownames(scz_meta))]
    allmeta[,3] = bd_meta$beta[match(all_genes, rownames(bd_meta))]
    allmeta[,4] = mdd_meta$beta[match(all_genes, rownames(mdd_meta))]
    allmeta[,5] = aad_meta$beta[match(all_genes, rownames(aad_meta))]
    allmeta[,6] = ibd_meta$beta[match(all_genes, rownames(ibd_meta))]
    allmeta[,7] = PTSD_meta$beta[match(all_genes, PTSD_meta$Geneid)]
    allmeta[,8] = PTSD_meta$beta[match(all_genes, MDD_meta$Geneid)]

    colnames(allmeta) = c("ASD", "SCZ", "BD", "MDD", "AAD", "IBD", paste0('PTSD/', area), paste0('MDD/', area))
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
    barplot$p.fdr = p.adjust(barplot$p.fdr,method="fdr")
    barplot$p.bootstrap = null$prob[findInterval(barplot$Mean, null$null)]
    barplot$p.symbol = ""
    barplot$p.symbol[barplot$p.bootstrap<0.05] = "*"
    barplot$p.symbol[barplot$p.bootstrap<0.01] = "**"
    barplot$p.symbol[barplot$p.bootstrap<0.001] = "***"
    barplot$Comparison = rownames(barplot)
    barplot$modality="microarray"
    barplot <- barplot[grep('(PTSD|MDD)/', barplot$Comparison), ]
    
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
}    
pdf(paste0(Figures_path, '/compare_with_gandall_bar_graph_with_MDD', '.pdf'))
plots
dev.off()

    


