###### --------------------------------------------------------------
# Library import
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(ggpubr)
library(data.table)

## format data with 2 sig figs
scaleFUN <- function(x) sprintf("%.1f", x)

## sort by GCS average strength, add column of index
GCS_strengthSort <- function(df){
  df%>%
    select('pos','strength_rep1','strength_rep2','strength_rep3')%>%
    filter(strength_rep1 < 100 & strength_rep2 < 100 & strength_rep3 < 100)%>%
    mutate(strength_mean = as.numeric((strength_rep1 +strength_rep2 + strength_rep3)/3))%>%
    arrange(strength_mean) %>%
    mutate(idx =  1:n())
}

Strength_ranked_plot <- function(sorted_df, figname){
  pt <- ggplot(sorted_df)+
    geom_point(aes(x=idx,y=strength_mean))+
    theme_bw()+
    ylab('Genomic Cleavage Strength')+
    xlab('')+
    geom_hline(yintercept=1, linetype="dashed", color = "blue")+
    geom_text(aes(100,1,label = 'y=1', vjust = -1),size=4)+
    ylim(0,45)+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.6))+
    ggtitle(figname)
  print(pt)
  return(pt)
}

Strength_ranked_plot_log <- function(sorted_df, figname){
  pt <- ggplot(sorted_df)+
    geom_point(aes(x=idx,y=strength_mean))+
    scale_y_log10(breaks = c(0,1,5,10,15,40))+
    theme_bw()+
    xlab('')+
    ylab('')+
    geom_hline(yintercept=1, linetype="dashed", color = "blue")+
    ylab(bquote(atop('Genomic Cleavage Strength', '(log scale)')))+
    theme(axis.text.x = element_text(angle = 0, hjust = 0.7),axis.title.y = element_text(size=10))
  print(pt)
  return(pt)
}


###### input --------------------------------------------------------------
dataDir <- '/GCS_calling/GCS_extended_version/'
dfFiles <-  list.files(file.path(getwd(), dataDir), pattern="*.txt", full.name=T)
names(dfFiles) <- gsub(".txt","",basename(dfFiles))
df_list <- lapply(dfFiles, read.table, sep = ',',header = TRUE) #list of GCSs from 5 FQs


# ranked order plot -------------------------------------------------------
df_sort = lapply(df_list, GCS_strengthSort) # sort the df lists


pt_1 = Strength_ranked_plot(df_sort[[1]],names(df_sort)[1])
pt_2 = Strength_ranked_plot(df_sort[[2]],names(df_sort)[2])
pt_3 = Strength_ranked_plot(df_sort[[3]],names(df_sort)[3])
pt_4 = Strength_ranked_plot(df_sort[[4]],names(df_sort)[4])
pt_5 = Strength_ranked_plot(df_sort[[5]],names(df_sort)[5])

lm_sum = ggarrange(pt_4,pt_3,pt_2,pt_1,pt_5, nrow =3, ncol=2, common.legend = TRUE, legend="bottom")
ggsave("GO_cutoff.png", units="in", lm_sum, width=7.5, height=9, dpi=1200)


for (i in 1:5){
  pt2 = Strength_ranked_plot_log(df_sort[[i]],names(df_sort)[i])
  ggsave(paste0(names(df_sort)[i],"_logscale.png"), units="in", pt2, width=2.9, height=2.2, dpi=1200)
}

