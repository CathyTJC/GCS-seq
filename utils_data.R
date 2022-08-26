library(DESeq2)
library(tidyverse)
library(RColorBrewer)
library(ggplot2)
library(ggrepel)
library(rtracklayer)
library(wesanderson)


## Read in mu-origin 1 annotation file, select needed columns and ouput a dataframe
readAnnotaion <- function(AnnotationPath =.anoPath,col_to_select=.colselect ){
  GTFannot<- rtracklayer::import(file.path(getwd(),AnnotationPath))
  gtf_df=as.data.frame(GTFannot) %>% ##a fuller version of the gtf
    select(col_to_select)%>% 
    unique()
  return(gtf_df)
}

## read in GCS and cleavage strength data 
# note 1 GCS has outlier strength already removed
readGCSStrength <- function(GCSPath = .gcspath, drug,suffix = .suffix,move_first=TRUE,col_to_select=FALSE){
  file_path <- file.path(getwd(),GCSPath,paste0(drug,suffix))
  if (move_first==TRUE){
    df <- read.table(file_path,sep = ',',header = TRUE)%>%
      select(-X) #remove first column (no info in the first column)
  }
  else{
    df <- read.table(file_path,sep = ',',header = TRUE)
  }
  if (col_to_select !=FALSE){
    df <-df[,col_to_select]
  }
  return(df)
}

## get counts from RNA-seq for each time point and drug for each replicates(3 total)
getRNAseqcounts <- function(RNAseqPath=.RNAseqpath,condition,drug){
  file_path <- file.path(getwd(),RNAseqPath,condition,drug)
  dfFiles =list.files(file_path, pattern = '.txt',full.names =TRUE)
  if (!is.null(dfFiles)){
    ldf <- lapply(dfFiles, read_table) }
  df<-left_join(ldf[[1]], ldf[[2]], by = 'Geneid')%>% 
    left_join(., ldf[[3]], by = 'Geneid')
  colnames(df) <- c('gene_ID',paste0(condition,'_',drug,'_1'),paste0(condition,'_',drug,'_2'),paste0(condition,'_',drug,'_3'))
  return(df)
}

## Add annotation according to position/coordinates
###### add annotation to ORF
#input: 1. GCS dataframe (first column as position) 2. annotation file
add_ano<-function(GCS,muori_ref){
  gene_name = c()
  gcs_set = c()
  for (site in GCS$pos){
    max_start_inx = max(which((muori_ref$start<=site)))
    if (muori_ref$end[max_start_inx]>=site){
      gene_name= c(gene_name, muori_ref$gene_name[max_start_inx])
      gcs_set= c(gcs_set, site)
    }
  }
  df<-data.frame(pos = gcs_set, gene_name=gene_name) %>%
    left_join(GCS, by = 'pos')%>%
    left_join(muori_ref, by = 'gene_name')%>%
    mutate(gene_id2 =gsub("^.*-", "", gene_id ))
  return(df)
}

#PCA
.PCA<-function(rld){
  p<-plotPCA(rld, intgroup="sampletype")
  print(p)
  p2<-plotPCA(rld, intgroup=c("sampletype","replicate"))
  print(p2)
  
  ###  PC3 and PC4 values for input to ggplot
  rld_mat <- assay(rld)
  pca <- prcomp(t(rld_mat))
  df <- cbind(meta, pca$x)
  ggplot(df) + 
    geom_point(aes(x=PC3, y=PC4, color = sampletype))
}

#
.cor_heatmap <- function(rld){
  rld_mat <- assay(rld)
  rld_cor <- cor(rld_mat)    ## cor() is a base R function
  heat.colors <- RColorBrewer::brewer.pal(6, "Blues")
  pheatmap(rld_cor, annotation = meta, color = heat.colors, border_color=NA, fontsize = 10, 
           fontsize_row = 10, height=15)
}


# combined_heatmap_view ---------------------------------------------------
#2 time points(2 combined df input), 3 replicate, a
.concat <- function(df_treat=df_treatment,df_recov=df_recovery,split = '2'){
  df_treat_combined <-left_join(df_treat[[1]], df_treat[[2]], by = 'gene_ID') %>% 
    left_join(., df_treat[[3]], by = 'gene_ID')%>% 
    left_join(., df_treat[[4]], by = 'gene_ID')
  df_recov_combined <-left_join(df_recov[[1]], df_recov[[2]], by = 'gene_ID') %>% 
    left_join(., df_recov[[3]], by = 'gene_ID')%>% 
    left_join(., df_recov[[4]], by = 'gene_ID') 
  
  if (split =='2'){
  dfcombined <- left_join(df_treat_combined, df_recov_combined,by='gene_ID') %>% 
    column_to_rownames(var="gene_ID") %>%
    as.matrix()}
  
  if (split =='recov'){
    dfcombined <- df_recov_combined %>% 
      column_to_rownames(var="gene_ID") %>%
      as.matrix()}
  
  if (split =='treat'){
    dfcombined <- df_treat_combined %>% 
      column_to_rownames(var="gene_ID") %>%
      as.matrix()}
  
  return(dfcombined)
}



## read all txt files at a directory
readfile<- function(inputPath = .path, pattern = '.txt',row_names = TRUE){
  file_path <- file.path(getwd(),inputPath )
  dfFiles =list.files(file_path, pattern = pattern,full.names =TRUE)
  if (!is.null(dfFiles)){
    if (row_names == TRUE){
      ldf <- lapply(dfFiles, read.table,row.names = 1) }
    else{ldf <- lapply(dfFiles, read.table)}
    full_name <- list.files(file_path, pattern = pattern)
    ls_name <- gsub(pattern, '',full_name )
    names(ldf) <- ls_name
    }
  return(ldf)
}

## select significant genes 
#(df$padj < 0.05)& (abs(df$log2FoldChange) >= 0.58),
.getSigGene<- function(df, .log2FoldChange=0.58, .padj=0.05, gene = gene){
  sigGene <- df[(df$padj < .padj)& (abs(df$log2FoldChange) >= .log2FoldChange),]$gene
  return(sigGene)
}

## select sig genes from multiple dataframes, output gene list
getsigGenes<-function(ldf){
  sigGene_list <- lapply(ldf, .getSigGene)
  return(sigGene_list)
}

## unravel the gene list to 1 single list
# read in the files, extract sig genes from each list and unlist
unwrap_genelist <- function(inputPath,log2FoldChange=0.58, padj=0.05, gene = gene,pattern=.pattern){
  ldf<-readfile(inputPath,pattern)
  sigGene_list <- getsigGenes(ldf)
  sigGenes <- unique(unlist(sigGene_list))
  return(sigGenes)
}

#prepare the df for heatmap
prepare_rnaseq_df <- function(df_treat=df_treatment,df_recov=df_recovery,inputPath,split = '2',log2FoldChange=0.58, padj=0.05, 
                       gene = gene, pattern=.pattern){
  df_combined <- .concat(df_treat,df_recov)
  sigGenes = unwrap_genelist(inputPath,log2FoldChange, padj, gene,pattern) #total of 3668 genes are dif expressed out of 4420
  sigdf <-df_combined[rownames(df_combined) %in% sigGenes,]
  return(sigdf)
}






