##############
#Gene enrichment plot
#######
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(dplyr)
library(ggpmisc)
library(ggpubr)

enrich_plot<-function(LEVO_treat,cutoff,figname){
  ## split GO term of GO number and detailed annotation
  LEVO_treat$name <-gsub("^.*~", "", LEVO_treat$Term)
  ## order the df by Fold.Enrichemnt
  LEVO_treat_or <- LEVO_treat[order(LEVO_treat$Fold.Enrichment,decreasing = TRUE),]%>% 
    as_tibble()
  LEVO_treat_or <- LEVO_treat_or[LEVO_treat_or$FDR<0.05,]
  
  m<- ggplot(LEVO_treat_or,aes(y= -log10(FDR), x = reorder(name, -log10(FDR)),color = Count, size=Count))+
    geom_point()+
    scale_color_continuous(low="blue", high="red",  limits=c(0, 500), breaks=seq(100, 400, by=100))+
    theme_bw()+
    coord_flip()+
    xlab ('')+
    ylim(0,17.5)+
    ggtitle(figname)+
    ylab(bquote(-log[10]~"FDR"))+
    scale_size_continuous(limits=c(0, 500), breaks=seq(100, 400, by=100))+
    guides(color= guide_legend(), size=guide_legend())
  return(m)
}

# GCS strength > 1 ----------------------------------------------------------------
MOXI <-read.table(file = paste0(getwd(), "/GO_ontology_data/GCS_s1/MOXI1_ano_chart.txt"),sep = '\t',header = TRUE)
Mplot <- enrich_plot(MOXI,0.05,'MOXI')

LEVO <-read.table(file = paste0(getwd(), "/GO_ontology_data/GCS_s1/LEVO1_ano_chart.txt"),sep = '\t',header = TRUE)
LEVOplot <- enrich_plot(LEVO,0.05,'LEVO')

GEMI <-read.table(file = paste0(getwd(), "/GO_ontology_data/GCS_s1/GEMI1_ano_chart.txt"),sep = '\t',header = TRUE)
Gplot <- enrich_plot(GEMI,0.05,'GEMI')

CIP <-read.table(file = paste0(getwd(), "/GO_ontology_data/GCS_s1/CIP1_ano_chart.txt"),sep = '\t',header = TRUE)
Cplot <- enrich_plot(CIP,0.05,'CIP')

NOR<-read.table(file = paste0(getwd(), "/GO_ontology_data/GCS_s1/NOR1_ano_chart.txt"),sep = '\t',header = TRUE)
Nplot <- enrich_plot(NOR,0.05, 'NOR')

lm_sum = ggarrange(Mplot,LEVOplot,Gplot,Cplot, nrow =2, ncol=2, common.legend = TRUE, legend="bottom")
ggsave("GO_enrichment_s1.pdf", units="in", lm_sum, width=7, height=7.5, dpi=1200)

# ALL_GCSs ----------------------------------------------------------------
enrich_plot_all<-function(LEVO_treat,cutoff,figname){
  ## split GO term of GO number and detailed annotation
  LEVO_treat$name <-gsub("^.*~", "", LEVO_treat$Term)
  ## order the df by Fold.Enrichemnt
  LEVO_treat_or <- LEVO_treat[order(LEVO_treat$Fold.Enrichment,decreasing = TRUE),]%>% 
    as_tibble()
  LEVO_treat_or <- LEVO_treat_or[LEVO_treat_or$FDR<0.05,]
  
  m<- ggplot(LEVO_treat_or,aes(y= -log10(FDR), x = reorder(name, -log10(FDR)),color = Count, size=Count))+
    geom_point()+
    scale_color_continuous(low="blue", high="red",  limits=c(0, 900), breaks=seq(100, 900, by=200))+
    theme_bw()+
    coord_flip()+
    xlab ('')+
    ylim(0,17.5)+
    ggtitle(figname)+
    ylab(bquote(-log[10]~"FDR"))+
    scale_size_continuous(limits=c(0, 900), breaks=seq(100, 900, by=200))+
    guides(color= guide_legend(), size=guide_legend())
  return(m)
}

MOXI <-read.table(file = paste0(getwd(), "/GO_ontology_data/allGCS/MOXI_ano_chart.txt"),sep = '\t',header = TRUE)
Mplot <- enrich_plot_all(MOXI,0.05,'MOXI')

LEVO <-read.table(file = paste0(getwd(), "/GO_ontology_data/allGCS/LEVO_ano_chart.txt"),sep = '\t',header = TRUE)
LEVOplot <- enrich_plot_all(LEVO,0.05,'LEVO')

GEMI <-read.table(file = paste0(getwd(), "/GO_ontology_data/allGCS/GEMI_ano_chart.txt"),sep = '\t',header = TRUE)
Gplot <- enrich_plot_all(GEMI,0.05,'GEMI')

CIP <-read.table(file = paste0(getwd(), "/GO_ontology_data/allGCS/CIP_ano_chart.txt"),sep = '\t',header = TRUE)
Cplot <- enrich_plot_all(CIP,0.05,'CIP')

NOR<-read.table(file = paste0(getwd(), "/GO_ontology_data/allGCS/NOR_ano_chart.txt"),sep = '\t',header = TRUE)
Nplot <-enrich_plot_all(NOR,0.05, 'NOR')

lm_sum = ggarrange(Mplot,LEVOplot,Gplot,Cplot, nrow =2, ncol=2, common.legend = TRUE, legend="bottom")
ggsave("GO_enrichment.pdf", units="in", lm_sum, width=9, height=7.5, dpi=1200)
