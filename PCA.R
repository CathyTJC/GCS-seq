################
#PCA analysis
#author: "Juechun Cathy Tang"
#date: "May-2022" updated

###### 
# Data Import
library(tidyverse)
library(reshape2)
library(ggplot2)
library(corrr)
library(gridExtra)
library(RColorBrewer)
library(gtools)


precheck<-function (trans_enrich){
  # Plot covariation between between samples
  cov<-trans_enrich %>% 
    ggplot(aes(LEVO_rep1, LEVO_rep2)) + geom_point() +
    geom_abline(colour = "brown")+
    scale_x_log10()+
    scale_y_log10()
  print(cov)
  
  #######
  #Plot Correlations between replicates
  trans_cts_corr <- trans_enrich %>% 
    # remove the column "pos", which we do not want to calculate correlation on
    select(-pos) %>% 
    # we use Spearman's correlation, a non-parametric metric based on ranks
    cor(method = "spearman")
  # Visualise the correlations between the first 5 samples
  coor<-rplot(trans_cts_corr) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(coor)

  #######
  ## log transformation of our cleavage strength data for clustering 
  # and visulization
  log_enrich <- trans_enrich %>% 
    # remove the column "gene"
    select(-pos) %>% 
    +1%>% 
    log(2) %>%
    mutate(pos = trans_enrich$pos)
  
  ## With the log transormed data, let's redo the above visualization
  # "gather" the counts data
  trans_cts_long <- log_enrich %>% 
    pivot_longer(cols = MOXI_rep1:NOR_rep3, 
                 names_to = "sample", 
                 values_to = "cts")
  
  enrich_dis<-trans_cts_long %>%
    ggplot(aes(cts, colour = sample)) + 
    theme_bw()+
    xlab('Cleavage Strength (log-tansformed)')+
    geom_freqpoly(binwidth =0.5) 
  plot(enrich_dis)
  
  return(log_enrich)
}

PCA<-function (df){
  pca_matrix <- df %>% 
    remove_rownames()%>% 
    # make the "track （GCS coordinate" column become the rownames of the table
    column_to_rownames("pos") %>% 
    # coerce to a matrix
    as.matrix() %>% 
    # transpose the matrix so that rows = samples and columns = variables
    t()
  
  # Perform the PCA, scale by unit variance
  pca <- prcomp(pca_matrix,scale. = TRUE)
  summary(pca)
  
  pc_eigenvalues <- pca$sdev^2
  pc_eigenvalues <- tibble(PC = factor(1:length(pc_eigenvalues)), 
                           variance = pc_eigenvalues) %>% 
    # add a new column with the percent variance
    mutate(pct = variance/sum(variance)*100) %>% 
    # add another column with the cumulative variance explained
    mutate(pct_cum = cumsum(pct))
  
  var_pc<-pc_eigenvalues %>% 
    ggplot(aes(x = PC)) +
    geom_col(aes(y = pct)) +
    geom_line(aes(y = pct_cum, group = 1)) + 
    geom_point(aes(y = pct_cum)) +
    labs(x = "Principal component", y = "Fraction variance explained")+
    theme_classic()
  print(var_pc)
  

  #PC scores
  pc_bplot<-pca$x %>% 
    # convert it to a tibble
    as_tibble(rownames = "sample") %>% 
    separate(col = sample,into = c("Treatment", "rep"), sep = "_") 
  
  pc_bplot$Drug<- factor(pc_bplot$Treatment,levels = c('MOXI','LEVO','GEMI','CIP','NOR'))
  
  # make the plot
  pc_bplot1 <- pc_bplot %>%
    ggplot(aes(x = PC1, y = PC2, color = Drug)) +
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)+
    theme_classic()+
    geom_point(size=3)
  plot(pc_bplot1)
  

  pc_bplot2 <- pc_bplot %>%
    ggplot(aes(x = PC1, y = PC3, color = Drug)) +
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)+
    theme_classic()+
    geom_point(size=3)
  plot(pc_bplot2)
  
  pc_bplot3 <- pc_bplot %>%
    ggplot(aes(x = PC2, y = PC3, color = Drug)) +
    scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)+
    theme_classic()+
    geom_point(size=3)
  plot(pc_bplot3)
  
  return(pca)
}

PCA_analysis<- function(pca){
  # Additional analysis: which genes have the most influence on each PC axis?
  pc_loadings <- pca$rotation %>% 
    as_tibble(rownames = "track")
  
  # There are over 7000 genes, select top genes and take a look
  top_GCS <- pc_loadings %>% 
    # select only the PCs we are interested in
    select(track, PC1, PC2)%>%
    # convert to a "long" format
    pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
    # for each PC
    group_by(PC) %>% 
    # arrange by descending order of loading
    arrange(desc(abs(loading))) %>% 
    # take the 10 top rows
    dplyr::slice(1:10) %>% 
    # pull the gene column as a vector
    pull(track) %>% 
    # ensure only unique genes are retained
    unique()
  
  top_loadings <- pc_loadings %>% 
    filter(track %in% top_GCS)
  
  LD_plot <- ggplot(data = top_loadings) +
    geom_segment(aes(x = 0, y = 0, xend = PC1, yend = PC2), 
                 arrow = arrow(length = unit(0.1, "in")),
                 colour = "brown") +
    geom_text(aes(x = PC1, y = PC2, label = track),
              position=position_jitter(width=0.005,height=0.005), size = 3) +
    scale_x_continuous(expand = c(0.02, 0.02))+
    ggtitle('Top GCS in PC1(x) and 2 (y)')
  print(LD_plot)
  
  ######
  # histogram of PC loadings
  # loading 1
  top_GCS1 <- pc_loadings %>% 
    # select only the PCs we are interested in
    select(track, PC1)%>%
    # convert to a "long" format
    pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
    # for each PC
    group_by(PC) %>% 
    # arrange by descending order of loading
    arrange(desc(abs(loading)))
  
  PC1_his<-ggplot(top_GCS1)+
    geom_histogram(mapping = aes(loading),bins = 100)
  print(PC1_his)

  #ggplot(top_GCS1)+
    #geom_histogram(mapping = aes(loading),bins = 100)
  
  print(top_GCS1$track[1:10])
  
  # loading 2
  top_GCS2 <- pc_loadings %>% 
    select(track, PC2)%>%
    pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
    group_by(PC) %>% 
    arrange(desc(abs(loading)))
  
  PC2_his<-ggplot(top_GCS2)+
    geom_histogram(mapping = aes(abs(loading)),bins = 100)
  print(PC2_his)
  
  # loading 3
  top_GCS3 <- pc_loadings %>% 
    select(track, PC3)%>%
    pivot_longer(matches("PC"), names_to = "PC", values_to = "loading") %>% 
    group_by(PC) %>% 
    arrange(desc(abs(loading)))
  PC3_his<-ggplot(top_GCS3)+
    geom_histogram(mapping = aes(abs(loading)),bins = 100)
  print(PC3_his)
  
  return(pc_loadings)
}


# Data input, PCA --------------------------------------------------------------
# Read in merged csv containing the cleavage strength for each replicate across 5 FQ treatment 


dataDir ='data/GCS_calling/GCS_merge_pseudo_strength.txt' # dataframe contains pseudo strengths, takes the union of all GCS
GCS_union_rep<-read.csv(file = file.path(getwd(), dataDir),header = TRUE)
pca<-PCA(GCS_union_rep)


# dataDir ='data/GCS/PCA/GCS_merge_shared.csv'
# GCS_shared <-read.csv(file = file.path(getwd(), dataDir),header = TRUE,row.names = 1) 
# #PCA of shared GCS
# pca<-PCA(GCS_shared)




