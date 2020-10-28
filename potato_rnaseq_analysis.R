library(tidyverse)
library(DESeq2)
library(RColorBrewer)
library(pheatmap)
library(goseq)
library(biomaRt)
library(gplots)
library(UpSetR)
library(EnhancedVolcano)

###########################################################################
#          Create output directories                                     #
###########################################################################
outdir <- "RESULTS"
dir.create(outdir)
dir.create(paste0(outdir,"PVALUE_HIST"))
#dir.create(paste0(outdir,"BOXPLOTS"))
dir.create(paste0(outdir,"DE_TABLES_SIGNIFICANT"))
dir.create(paste0(outdir,"DE_TABLES_FULL"))
dir.create(paste0(outdir,"EXPLORATORY"))
dir.create(paste0(outdir,"OBJECTS"))
dir.create(paste0(outdir,"GO_TABLES_SIGNIFICANT"))
dir.create(paste0(outdir,"GO_TABLES_FULL"))
dir.create(paste0(outdir,"HEATMAPS"))
dir.create(paste0(outdir,"VOLCANO_PLOTS"))

#################################################################################
# Read input files                                                              #
#################################################################################
counts_df <- read.csv("all_comparisons_merged.csv")

head(counts_df)
dim(counts_df)

summary(counts_df$mapping_type)
#cluster_bed     genomic 
#2171       27398 

#select raw counts
raw_counts <- counts_df[,c(1,28:39)]
dim(raw_counts)
#29569 12

# set gene id as rownames
rownames(raw_counts) <- raw_counts$gene_id
raw_counts <- raw_counts[,-1]
head(raw_counts)


# sample, condition
coldata <- read.csv("coldata.csv")
coldata2 <- coldata %>% dplyr::select(sample, condition)
dim(coldata2)
#[1] 12  2

rownames(coldata2) <- coldata2$sample
coldata2<- data.frame(coldata2)
coldata2
dim(coldata2)
#[1] 12 2


######################################################################
#     Create deseq2 object from input files and prefilter            #  
#######################################################################
####create deseq2 object
dds <- DESeqDataSetFromMatrix(countData = raw_counts,
                              colData = coldata2,
                              design= ~ condition)

nrow(dds)
#remove rows of the DESeqDataSet that have no counts, or only a single count across all samples.
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)
#[1] 27426 

#######################################################################
#        Exploratory analysis                                         #
#######################################################################
# do pca 
vsd <- vst(dds, blind=FALSE)
#pca plot
tiff(paste0(outdir,"EXPLORATORY/","pca.tiff"),height = 5, width = 5, units = "in", res = 300)
plotPCA(vsd, intgroup="condition")
dev.off()

##heatmap of sample dists
sampleDists <- dist(t(assay(vsd)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(vsd$condition, vsd$type, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)

tiff(paste0(outdir,"EXPLORATORY/","Sampledists.tiff"),height = 5, width = 5, units = "in", res = 300)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors)
dev.off()

#mds plot
tiff(paste0(outdir,"EXPLORATORY/","mdsplot.tiff"),height = 5, width = 5, units = "in", res = 300)
mds <- as.data.frame(colData(vsd))  %>%
  cbind(cmdscale(sampleDistMatrix))
ggplot(mds, aes(x = `1`, y = `2`, color = condition)) +
  geom_point(size = 3) + coord_fixed() + ggtitle("MDS with VST data")
dev.off()

#######################################################################
# Run DeSeq2                                                          #
#######################################################################
#run deseq2
dds <- DESeq(dds)
saveRDS(dds, paste0(outdir,"OBJECTS/","deseq-dds.rds"))

resultsNames(dds) # lists the coefficients


###############################################################################
#                  Get gene annotation and GO annotation                      #
###############################################################################
#get plant ensembl info
listMarts(host="plants.ensembl.org")
db = useMart("plants_mart", host="plants.ensembl.org")
db

#see all datasets in this mart
datasets <- listDatasets(db)
datasets

#use the potato genome dataset
#stuberosum_eg_gene
# stuberosum_eg_gene     Solanum tuberosum genes (SolTub_3.0)
ensembl = useDataset("stuberosum_eg_gene",mart= db)


#get some features 
features = getBM(attributes = c("ensembl_gene_id",
                                "description","chromosome_name"
                                ,"start_position", "end_position"),
                 filters = c("ensembl_gene_id"),
                 values = rownames(dds), mart = ensembl)
dim(features)
#[1] 24593 5 

#add gene length info to features df
features$gene_length <- features$end_position - features$start_position
#write to csv
write.csv(features, paste0(outdir,"/features.csv"))

#######################################################
#see all possible features that I can get 
listAttributes(ensembl)

##get GO information for ensembl genes
features_2 <- getBM(attributes=c('ensembl_gene_id', 'description', 'chromosome_name', 'start_position', 'end_position', 'kegg_enzyme','plant_reactome_pathway', 'go_id', 'go_linkage_type', 'name_1006', 'namespace_1003','definition_1006'), values = rownames(dds), mart = ensembl)
head(features_2)
dim(features_2)
#[1] 131340 11  #genes listed multiple times if have many go terms 
write.csv(features_2, "features_2.csv")

unique(features_2$go_linkage_type)
#[1] ""    "IEA" "IDA" "IPI" "IMP" "ISS" "IC"  "TAS" "RCA" "IEP" "NAS" "ND" 

unique(features_2$namespace_1003)
#[1] ""                   "biological_process" "molecular_function" "cellular_component" "go" 

#######
features_3 <- getBM(attributes=c('ensembl_gene_id', 'description'), values = rownames(dds), mart = ensembl)
head(features_3)


######################################################################################
#                       Report Results                                               # 
######################################################################################

##function to generate tables of deseq2 results, produce histogram of pvalues
deseq_results <- function(dds, case, control, alpha=0.05, features_df){
  
  require(dplyr) 
  require(tidyverse)
  res <- results(dds, contrast = c("condition", case, control), alpha = alpha)
  outputname <- paste0(case, "_vs_", control)
  
  #add gene description 
  res_df<- data.frame(res)
  res_df$gene_id <- rownames(res_df)
  res_df_2 <- left_join(res_df, features_df, by = c("gene_id" = "ensembl_gene_id"))
  
  #sbset to significant results
  sigres <- res_df_2[which(res_df_2$padj < alpha), ]
  print(dim(sigres))
  print(summary(sigres$log2FoldChange))
  
  #keep only significant results with absolute log2 fold change >2 
  sigres_2 <- sigres[which(abs(sigres$log2FoldChange) >= 2), ]
  print(summary(sigres_2$log2FoldChange))
  print(dim(sigres_2))
  
  #histogram of pvalues 
  pdf(paste0(outdir,"PVALUE_HIST/",outputname,"-pvalue-hist.pdf"))
  hist(res_df$pvalue[res$baseMean > 1], col = "darkred", main = outputname, xlab = "p-values")
  dev.off()
  
  #put back in rownames
  rownames(sigres)<-sigres$gene_id
  rownames(sigres_2)<-sigres_2$gene_id
  rownames(res_df_2)<-res_df_2$gene_id
  
  #write significant results (padjust < alpha)
  write.csv(sigres, paste0(outdir,"DE_TABLES_SIGNIFICANT/",outputname,"-sigres.csv"))
  #wite significant results with padjust < alpha and absolute log2 fold change >2 
  write.csv(sigres_2, paste0(outdir,"DE_TABLES_SIGNIFICANT/",outputname,"-sigres_log2FC_filt.csv"))
  #all results
  write.csv(res_df_2, paste0(outdir,"DE_TABLES_FULL/",outputname,"-allres.csv"))
  
  #return both
  df_list <- list(sig_df = sigres, sig_df_logfilt = sigres_2)
  
  
}


###run analysis   (# FDR < 0.05, log2FC >= 2 or log2FC <= 2)
#2 dataframes stored in a list for each result variable below
###
l_vs_c<- deseq_results(dds = dds, case = "L", control = "C", alpha = 0.05, features_df = features_3) #?
f_vs_c <- deseq_results(dds = dds, case = "F", control = "C", alpha = 0.05, features_df = features_3) #279
m_vs_c<- deseq_results(dds = dds, case = "M", control = "C", alpha = 0.05, features_df = features_3) #19

##
l_vs_f <-deseq_results(dds = dds, case = "L", control = "F", alpha = 0.05, features_df = features_3) #3
l_vs_m <-deseq_results(dds = dds, case = "L", control = "M", alpha = 0.05, features_df = features_3) #243
f_vs_m <- deseq_results(dds = dds, case = "F", control = "M", alpha = 0.05, features_df = features_3) #48 

#access elements in list example:
#l_vs_c$sig_df
#l_vs_c$sig_df_logfilt


###################################################################################
#                      Volcano Plots                                              #
###################################################################################


##number of up and down de genes in each condition 
#triangle type plot

tiff(paste0(outdir,"VOLCANO_PLOTS/L_vs_C-allde.tiff"),height = 7, width = 7, units = "in", res = 300)
EnhancedVolcano(l_vs_c$sig_df,
                lab = rownames(l_vs_c$sig_df),
                x = 'log2FoldChange',
                y = 'padj',
                title = "L vs C",
                subtitle = "",
                xlim = c(-8, 8), 
                ylab = "-Log10 Ajusted P", 
                pCutoff = 0.05, 
                FCcutoff = 2)
dev.off()


tiff(paste0(outdir,"VOLCANO_PLOTS/f_vs_C-allde.tiff"),height = 7, width = 7, units = "in", res = 300)
EnhancedVolcano(f_vs_c$sig_df,
                lab = rownames(f_vs_c$sig_df),
                x = 'log2FoldChange',
                y = 'padj',
                title = "F vs C",
                subtitle = "",
                ylab = "-Log10 Ajusted P", 
                pCutoff = 0.05, 
                FCcutoff = 2,
                xlim = c(-7, 7))
dev.off()



tiff(paste0(outdir,"VOLCANO_PLOTS/m_vs_C-allde.tiff"),height = 7, width = 7, units = "in", res = 300)
EnhancedVolcano(m_vs_c$sig_df,
                lab = rownames(m_vs_c$sig_df),
                x = 'log2FoldChange',
                y = 'padj',
                title = "M vs C",
                subtitle = "",
                ylab = "-Log10 Ajusted P", 
                pCutoff = 0.05, 
                FCcutoff = 2, 
                xlim = c(-5,5))
dev.off()

tiff(paste0(outdir,"VOLCANO_PLOTS/L_vs_F-allde.tiff"),height = 7, width = 7, units = "in", res = 300)
EnhancedVolcano(l_vs_f$sig_df,
                lab = rownames(l_vs_f$sig_df),
                x = 'log2FoldChange',
                y = 'padj',
                title = "L vs F",
                subtitle = "",
                ylab = "-Log10 Ajusted P", 
                pCutoff = 0.05, 
                FCcutoff = 2, 
                xlim = c(-4,4))
dev.off()


tiff(paste0(outdir,"VOLCANO_PLOTS/L_vs_M-allde.tiff"),height = 7, width = 7, units = "in", res = 300)
EnhancedVolcano(l_vs_m$sig_df,
                lab = rownames(l_vs_m$sig_df),
                x = 'log2FoldChange',
                y = 'padj',
                title = "L vs M",
                subtitle = "",
                ylab = "-Log10 Ajusted P", 
                pCutoff = 0.05, 
                FCcutoff = 2, 
                xlim = c(-6,6))
dev.off()

tiff(paste0(outdir,"VOLCANO_PLOTS/F_vs_M-allde.tiff"),height = 7, width = 7, units = "in", res = 300)
EnhancedVolcano(f_vs_m$sig_df,
                lab = rownames(f_vs_m$sig_df),
                x = 'log2FoldChange',
                y = 'padj',
                title = "F vs M",
                subtitle = "",
                ylab = "-Log10 Ajusted P", 
                pCutoff = 0.05, 
                FCcutoff = 2, 
                xlim =c (-6,6))
dev.off()

#################

###dispersion plot
tiff(paste0(outdir,"EXPLORATORY/dispresion-check.tiff"),height = 5, width = 5, units = "in", res = 300)
plotDispEsts(dds)
dev.off()


#write out normalised counts
norm_counts_df <- counts(dds, normalized = T)
write.csv(norm_counts_df, paste0(outdir,"/deseq2-normalised-counts.csv"))


###############################################################################
#                      GO Enrichment Analysis                                 #
###############################################################################

##make category mappings
go_category_info <- features_2 %>% dplyr::select(ensembl_gene_id, go_id)
kegg_category_info <- features_2 %>% dplyr::select(ensembl_gene_id, kegg_enzyme)
reactome_category_info <- features_2 %>% dplyr::select(ensembl_gene_id, plant_reactome_pathway)


run_go_enrichment<-function(dds, sig_res_df, name, features_df, features_df_2){
  
  #assign 1 to de genes and 0 to no de genes
  de_genes <- as.integer(rownames(dds) %in% rownames(sig_res_df))
  names(de_genes) <-rownames(dds)
  print(table(de_genes))
  
  
  #get only gene vector for names that have ensembl ids 
  de_genes_2 <- de_genes[names(de_genes) %in% features$ensembl_gene_id]
  length(de_genes_2)
  print(table(de_genes_2))
  
  #extract gene lengths 
  rownames(features) <- features$ensembl_gene_id
  #ensure order is the same 
  features_new <- features[names(de_genes_2),]
  
  #check order
  print(all(names(de_genes_2) == rownames(features_new)))
  stopifnot(all(names(de_genes_2) == rownames(features_new)))
  
  #get gene lengths
  de_genes_2_length <- features_new$gene_length
  names(de_genes_2_length) <- rownames(features_new)
  de_genes_2_length
  
  ###run goseq
  pwf = nullp(de_genes_2, bias.data = de_genes_2_length)
  go_wall_res <- goseq(pwf,gene2cat = go_category_info)
  
  #produced BH adjusted p value for overrep terms
  go_wall_res$over_padj <- p.adjust(go_wall_res$over_represented_pvalue, method = "BH")
  write.csv(go_wall_res, paste0(outdir,"GO_TABLES_FULL/", name, "-go_all.csv"))
  
  #get significant go terms
  go_wall_res_sig <- go_wall_res[go_wall_res$over_padj < 0.05,]
  
  #remove any empty category 
  go_wall_res_sig_2 <- go_wall_res_sig %>% filter(category != "")
  print(go_wall_res_sig_2)
  
  if(nrow(go_wall_res_sig_2) > 0){
    
    write.csv(go_wall_res_sig_2, paste0(outdir,"GO_TABLES_SIGNIFICANT/", name, "-go_sig-logfilt.csv"))
    #get all genes in these significant categories
    features_go_sub<- features_2[features_2$go_id %in% go_wall_res_sig_2$category,]
  
    #get only de gene names
    de_genes_only<- de_genes_2[de_genes_2==1]
    print(table(de_genes_only))
  
    #get all de_genes and associated go terms here 
    all_de_genes_with_go_terms<-features_go_sub[features_go_sub$ensembl_gene_id %in% names(de_genes_only),]
    write.csv(all_de_genes_with_go_terms, paste0(outdir,"GO_TABLES_SIGNIFICANT/", name, "-go_sig_de_gene_names-logfilt.csv"))
  }
}

#rn function
run_go_enrichment(dds=dds, sig_res_df = l_vs_c$sig_df_logfilt, "l_vs_c", features_df = features, features_df_2 = features_2)
run_go_enrichment(dds=dds, sig_res_df = f_vs_c$sig_df_logiflt, "f_vs_c", features_df = features, features_df_2 = features_2)                         
run_go_enrichment(dds=dds, sig_res_df = m_vs_c$sig_df_logfilt, "m_vs_c", features_df = features, features_df_2 = features_2) 

run_go_enrichment(dds=dds, sig_res_df = l_vs_f$sig_df_logfilt, "l_vs_f", features_df = features, features_df_2 = features_2) 
run_go_enrichment(dds=dds, sig_res_df = l_vs_m$sig_df_logfilt, "l_vs_m", features_df = features, features_df_2 = features_2) 
run_go_enrichment(dds=dds, sig_res_df = f_vs_m$sig_df_logfilt, "f_vs_m", features_df = features, features_df_2 = features_2) 




###############################################################################
#                      UpSetR Plots                                          #
###############################################################################
input_gene_list <- list(L_vs_C = l_vs_c$sig_df_logfilt$gene_id, 
                  F_vs_C = f_vs_c$sig_df_logfilt$gene_id,
                  M_vs_C = m_vs_c$sig_df_logfilt$gene_id,
                  L_vs_F = l_vs_f$sig_df_logfilt$gene_id, 
                  L_vs_M = l_vs_m$sig_df_logfilt$gene_id,
                  F_vs_M = f_vs_m$sig_df_logfilt$gene_id)


tiff(paste0(outdir,"UpsetR-de-genes-withemptyintersection-logfilt.tiff"), width = 11, height = 6, res = 300, units="in")
upset(fromList(input_gene_list), order.by = "freq", nsets = 6,  empty.intersections = "on", text.scale = 1.5)
dev.off()


tiff(paste0(outdir,"UpsetR-de-genes-noempty-logfilt.tiff"), width = 11, height = 6, res = 300, units="in")
upset(fromList(input_gene_list), order.by = "freq", nsets = 6, text.scale = 1.5)
dev.off()


#seperate up and down differentially expressed genes using log2FoldChange variable

#l vs c 
l_vs_c_up<- l_vs_c$sig_df_logfilt[l_vs_c$sig_df_logfilt$log2FoldChange >0,]
l_vs_c_down<- l_vs_c$sig_df_logfilt[l_vs_c$sig_df_logfilt$log2FoldChange <0,]

write.csv(l_vs_c_up, "l_vs_c-sigres_log2FC_filt_UP.csv" )
write.csv(l_vs_c_down, "l_vs_c_sigres_log2FC_filt_DOWN.csv")

nrow(l_vs_c$sig_df_logfilt)
nrow(l_vs_c_up) +nrow(l_vs_c_down)

#m vs c 
m_vs_c_up<- m_vs_c$sig_df_logfilt[m_vs_c$sig_df_logfilt$log2FoldChange >0,]
m_vs_c_down<- m_vs_c$sig_df_logfilt[m_vs_c$sig_df_logfilt$log2FoldChange <0,]

nrow(m_vs_c$sig_df_logfilt)
nrow(m_vs_c_up) +nrow(m_vs_c_down)

write.csv(m_vs_c_up, "m_vs_c-sigres_log2FC_filt_UP.csv")
write.csv(m_vs_c_down, "m_vs_c_sigres_log2FC_filt_DOWN.csv")


# f vs c 
f_vs_c_up<- f_vs_c$sig_df_logfilt[f_vs_c$sig_df_logfilt$log2FoldChange >0,]
f_vs_c_down<- f_vs_c$sig_df_logfilt[f_vs_c$sig_df_logfilt$log2FoldChange <0,]

nrow(f_vs_c$sig_df_logfilt)
nrow(f_vs_c_up) +nrow(f_vs_c_down)


write.csv(f_vs_c_up, "f_vs_c-sigres_log2FC_filt_UP.csv")
write.csv(f_vs_c_down, "f_vs_c_sigres_log2FC_filt_DOWN.csv")



# l vs f 

l_vs_f_up<- l_vs_f$sig_df_logfilt[l_vs_f$sig_df_logfilt$log2FoldChange >0,]
l_vs_f_down<- l_vs_f$sig_df_logfilt[l_vs_f$sig_df_logfilt$log2FoldChange <0,]

write.csv(l_vs_f_up, "l_vs_f-sigres_log2FC_filt_UP.csv")
write.csv(l_vs_f_down, "l_vs_f_sigres_log2FC_filt_DOWN.csv")

nrow(l_vs_f$sig_df_logfilt)
nrow(l_vs_f_up) +nrow(l_vs_f_down)

# L vs m 
l_vs_m_up<- l_vs_m$sig_df_logfilt[l_vs_m$sig_df_logfilt$log2FoldChange >0,]
l_vs_m_down<- l_vs_m$sig_df_logfilt[l_vs_m$sig_df_logfilt$log2FoldChange <0,]

nrow(l_vs_m$sig_df_logfilt)
nrow(l_vs_m_up) +nrow(l_vs_m_down)


write.csv(l_vs_m_up, "l_vs_m-sigres_log2FC_filt_UP.csv")
write.csv(l_vs_m_down, "l_vs_m_sigres_log2FC_filt_DOWN.csv")

# f vs m 
f_vs_m_up<- f_vs_m$sig_df_logfilt[f_vs_m$sig_df_logfilt$log2FoldChange >0,]
f_vs_m_down<- f_vs_m$sig_df_logfilt[f_vs_m$sig_df_logfilt$log2FoldChange <0,]

nrow(f_vs_m$sig_df_logfilt)
nrow(f_vs_m_up) +nrow(f_vs_m_down)


write.csv(f_vs_m_up, "f_vs_m-sigres_log2FC_filt_UP.csv")
write.csv(f_vs_m_down, "f_vs_m_sigres_log2FC_filt_DOWN.csv")


#combine into lists 

input_gene_list_up <- list(L_vs_C = l_vs_c_up$gene_id, 
                        F_vs_C = f_vs_c_up$gene_id,
                        M_vs_C = m_vs_c_up$gene_id,
                        L_vs_F = l_vs_f_up$gene_id, 
                        L_vs_M = l_vs_m_up$gene_id,
                        F_vs_M = f_vs_m_up$gene_id)



input_gene_list_down <- list(L_vs_C = l_vs_c_down$gene_id, 
                           F_vs_C = f_vs_c_down$gene_id,
                           M_vs_C = m_vs_c_down$gene_id,
                           L_vs_F = l_vs_f_down$gene_id, 
                           L_vs_M = l_vs_m_down$gene_id,
                           F_vs_M = f_vs_m_down$gene_id)


tiff(paste0(outdir,"UpsetR-de-genes-withemptyintersection-logfilt-UPREGULATED_GENES_ONLY.tiff"), width = 11, height = 6, res = 300, units="in")
upset(fromList(input_gene_list_up), order.by = "freq", nsets = 6,  empty.intersections = "on", text.scale = 1.5)
dev.off()


tiff(paste0(outdir,"UpsetR-de-genes-noempty-logfilt-UPREGULATED_GENES_ONLY.tiff"), width = 11, height = 6, res = 300, units="in")
upset(fromList(input_gene_list_up), order.by = "freq", nsets = 6, text.scale = 1.5)
dev.off()



tiff(paste0(outdir,"UpsetR-de-genes-withemptyintersection-logfilt-DOWNREGULATED_GENES_ONLY.tiff"), width = 11, height = 6, res = 300, units="in")
upset(fromList(input_gene_list_down), order.by = "freq", nsets = 6,  empty.intersections = "on", text.scale = 1.5)
dev.off()


tiff(paste0(outdir,"UpsetR-de-genes-noempty-logfilt-DOWNREGULATED_GENES_ONLY.tiff"), width = 11, height = 6, res = 300, units="in")
upset(fromList(input_gene_list_down), order.by = "freq", nsets = 6, text.scale = 1.5)
dev.off()



########################################################################################################
#                                     Heatmap of all DE genes                                          #
########################################################################################################
all_de_df <- rbind(l_vs_c$sig_df, f_vs_c$sig_df, m_vs_c$sig_df, l_vs_f$sig_df, l_vs_m$sig_df, f_vs_m$sig_df)
dim(all_de_df)
length(unique(all_de_df$gene_id))
#[1] 1939 unique de genes 

#check 
nrow(l_vs_c$sig_df) + nrow(f_vs_c$sig_df) + nrow(m_vs_c$sig_df)+nrow(l_vs_f$sig_df)+nrow(l_vs_m$sig_df)+nrow(f_vs_m$sig_df)


norm_de_df <- norm_counts_df[rownames(norm_counts_df) %in% all_de_df$gene_id, ]
dim(norm_de_df)
#[1] 1939 


#fix colnames 
head(norm_de_df)
colnames(norm_de_df)<- c("C1", "C2", "C3", "F1", "F2", "F3", "L1", "L2", "L3", "M1", "M2", "M3")

## Plot heatmap
tiff(paste0(outdir,"HEATMAPS/pheatmap-rowscale-pearsoncor-complete.tiff"),height = 5, width = 5, units = "in", res = 300)
pheatmap(norm_de_df, cluster_rows = T, show_rownames=F,scale = "row",clustering_distance_rows
= "correlation", clustering_distance_columns = "correlation", clustering_method = "complete",  col=rev(brewer.pal(11,"RdBu")))
dev.off()

########################################################################################################
#                                    Session Information                                               #
########################################################################################################
sessionInfo()

# R version 3.6.1 (2019-07-05)
# Platform: x86_64-w64-mingw32/x64 (64-bit)
# Running under: Windows 10 x64 (build 18362)
# 
# Matrix products: default
# 
# locale:
#   [1] LC_COLLATE=English_United Kingdom.1252  LC_CTYPE=English_United Kingdom.1252    LC_MONETARY=English_United Kingdom.1252
# [4] LC_NUMERIC=C                            LC_TIME=English_United Kingdom.1252    
# 
# attached base packages:
#   [1] parallel  stats4    stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] EnhancedVolcano_1.2.0       ggrepel_0.8.2               UpSetR_1.4.0                gplots_3.0.3               
# [5] biomaRt_2.40.5              goseq_1.36.0                geneLenDataBase_1.20.0      BiasedUrn_1.07             
# [9] pheatmap_1.0.12             RColorBrewer_1.1-2          DESeq2_1.24.0               SummarizedExperiment_1.14.1
# [13] DelayedArray_0.10.0         BiocParallel_1.18.1         matrixStats_0.55.0          Biobase_2.44.0             
# [17] GenomicRanges_1.36.1        GenomeInfoDb_1.20.0         IRanges_2.18.3              S4Vectors_0.22.1           
# [21] BiocGenerics_0.30.0         forcats_0.4.0               stringr_1.4.0               dplyr_0.8.3                
# [25] purrr_0.3.2                 readr_1.3.1                 tidyr_1.0.0                 tibble_2.1.3               
# [29] ggplot2_3.2.1               tidyverse_1.3.0            
# 
# loaded via a namespace (and not attached):
#   [1] colorspace_1.4-1         htmlTable_1.13.2         XVector_0.24.0           base64enc_0.1-3          fs_1.3.1                
# [6] rstudioapi_0.10          bit64_0.9-7              AnnotationDbi_1.46.1     lubridate_1.7.4          xml2_1.2.2              
# [11] splines_3.6.1            geneplotter_1.62.0       knitr_1.25               Formula_1.2-3            jsonlite_1.6            
# [16] Rsamtools_2.0.3          broom_0.5.2              annotate_1.62.0          GO.db_3.8.2              cluster_2.1.0           
# [21] dbplyr_1.4.2             compiler_3.6.1           httr_1.4.1               backports_1.1.5          assertthat_0.2.1        
# [26] Matrix_1.2-17            lazyeval_0.2.2           cli_1.1.0                acepack_1.4.1            htmltools_0.4.0         
# [31] prettyunits_1.0.2        tools_3.6.1              gtable_0.3.0             glue_1.3.1               GenomeInfoDbData_1.2.1  
# [36] Rcpp_1.0.2               cellranger_1.1.0         vctrs_0.3.0              Biostrings_2.52.0        gdata_2.18.0            
# [41] nlme_3.1-140             rtracklayer_1.44.4       xfun_0.10                rvest_0.3.5              lifecycle_0.1.0         
# [46] gtools_3.8.1             XML_3.98-1.20            zlibbioc_1.30.0          scales_1.0.0             hms_0.5.2               
# [51] memoise_1.1.0            gridExtra_2.3            rpart_4.1-15             latticeExtra_0.6-28      stringi_1.4.3           
# [56] RSQLite_2.1.2            genefilter_1.66.0        checkmate_1.9.4          caTools_1.17.1.2         GenomicFeatures_1.36.4  
# [61] rlang_0.4.5              pkgconfig_2.0.3          bitops_1.0-6             lattice_0.20-38          GenomicAlignments_1.20.1
# [66] htmlwidgets_1.5.1        bit_1.1-14               tidyselect_0.2.5         plyr_1.8.4               magrittr_1.5            
# [71] R6_2.4.0                 generics_0.0.2           Hmisc_4.3-1              DBI_1.0.0                mgcv_1.8-28             
# [76] pillar_1.4.2             haven_2.2.0              foreign_0.8-71           withr_2.1.2              survival_3.1-11         
# [81] RCurl_1.95-4.12          nnet_7.3-12              modelr_0.1.5             crayon_1.3.4             KernSmooth_2.23-15      
# [86] progress_1.2.2           locfit_1.5-9.1           grid_3.6.1               readxl_1.3.1             data.table_1.12.4       
# [91] blob_1.2.0               reprex_0.3.0             digest_0.6.21            xtable_1.8-4             munsell_0.5.0 


