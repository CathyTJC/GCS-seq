# This script performs differential gene expression analysis and saves the results
# Author: Cathy Tang

### Bioconductor and CRAN libraries used
source('utils_data.R')

# DEseq2 ------------------------------------------------------------------
# inputDF: gene_ID/treatment_negative_1/treatment_negative_2/treatment_negative_3
Deseq <-function(control, sample,GCS_LEVO, plot =FALSE,annotation=anoDF,drug_name){
  ### combine into count matrix_raw
  LEVO_neg_30min <-left_join(control, sample, by = 'gene_ID') %>% column_to_rownames(var="gene_ID") %>%
    as.matrix()
  ## Create a metadata
  sampletype <- factor(c(rep("untreated",3), rep("treated", 3)))
  replicate<- factor(rep(c('n1','n2','n3'),2))
  meta <- data.frame(sampletype, replicate,row.names = colnames(LEVO_neg_30min))
  ## Create DESeq2Dataset object
  dds <- DESeqDataSetFromMatrix(LEVO_neg_30min, colData = meta, design = ~ sampletype + replicate)
  dds <- estimateSizeFactors(dds)
  sizeFactors(dds)
  normalized_counts <- counts(dds, normalized=TRUE)
  #save normalized count table
  rld <- rlog(dds, blind=TRUE)
  
  ## Create DESeq2Dataset object， with betaPrior set already shunken
  dds <- DESeqDataSetFromMatrix(LEVO_neg_30min, colData = meta, design = ~ sampletype + replicate)
  ## Run analysis
  dds <- DESeq(dds, fitType = "parametric",minReplicatesForReplace = 7, betaPrior = TRUE)
  
  ## Define contrasts
  contrast_oe <- c("sampletype", "treated", "untreated")
  res <- results(dds, contrast=contrast_oe,independentFiltering =TRUE,cooksCutoff =TRUE)
  
  ## examine results
  ## Save the results to compare
  res_sorted <- res[order(res$padj), ]
  summary(res, alpha = 0.05)
  
  ## Obtain logical vector where TRUE values denote padj values < 0.05 and fold change > 1.5 in either direction
  res_tableOE_tb <- res_sorted %>% 
    data.frame() %>%
    rownames_to_column(var="gene") %>% 
    as_tibble() %>%
    left_join(annotation, by=c("gene" = "gene_id"))%>%
    mutate(threshold_OE = padj < 0.05 & abs(log2FoldChange) >= 0.58)
  #save sorted LFC results 
  write.table(res_tableOE_tb , file=file.path(getwd(),'data/gene_exp/differential_expression',paste0(drug_name, "_deseq.txt")), sep="\t", quote=F, col.names=NA)

  ##QC checks
  if (plot ==TRUE){
    ##Visualization
    .PCA(rld)
    .cor_heatmap(rld)
    plotDispEsts(dds)  # dispersion plots
    # MA plot using unshrunken fold changes
    plotMA(res,ylim = range(res$log2FoldChange, na.rm = TRUE))
  }
  return(list(res_tableOE_tb,meta,normalized_counts))
}




#read annoation file
.anoPath <- '/data/REFERENCE/GCF_000005845.2_ASM584v2_genomic_insert_hierexon.gffread.gtf'
.colselect = c("gene_id","gene_name",'width','start','end','strand')
anoDF<-readAnnotaion(AnnotationPath =.anoPath,col_to_select=.colselect)

#read GCS file
.gcspath <- '/data/GCS'
.suffix <- '_sh.txt'
MOXI <- readGCSStrength(GCSPath = .gcspath, drug = 'MOXI',suffix  =.suffix)
LEVO <- readGCSStrength(GCSPath = .gcspath, drug = 'LEVO',suffix  =.suffix)
MOXI_ano <-  add_ano(MOXI,anoDF) #add annotation to GCS coordinates
LEVO_ano <-  add_ano(LEVO,anoDF)

#read gene expression files (raw counts)
.RNAseqpath <- 'data/gene_exp/COUNTS/RAW'
.condition = c('treatment','recovery')
.drug = c('negative','MOXI','LEVO','MMC')
df_treatment <- lapply(.drug,getRNAseqcounts, RNAseqPath=.RNAseqpath, condition = .condition[1])
df_recovery <- lapply(.drug,getRNAseqcounts, RNAseqPath=.RNAseqpath, condition = .condition[2])


Treatment_LEVO <- Deseq(df_treatment[[1]],df_treatment[[3]],drug_name='Treatment_LEVO')
Treatment_MOXI <- Deseq(df_treatment[[1]],df_treatment[[2]],drug_name='Treatment_MOXI')
Treatment_MMC <- Deseq(df_treatment[[1]],df_treatment[[4]],drug_name='Treatment_MMC')
Recovery_LEVO <- Deseq(df_recovery[[1]],df_recovery[[3]],drug_name='Recovery_LEVO')
Recovery_MOXI <- Deseq(df_recovery[[1]],df_recovery[[2]],drug_name='Recovery_MOXI')
Recovery_MMC <- Deseq(df_recovery[[1]],df_recovery[[4]],drug_name='Recovery_MMC')




  
  
  
  


