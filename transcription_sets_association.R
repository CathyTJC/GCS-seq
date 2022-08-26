###############################################
##Juechun Tang, 2021##
##GCS-seq##
#Script visualizes GCSs within different sets of transcription level
###############################################

# Library Import -------------------------------------------------------------------
source('utils_data.R')
library(gridExtra)
library(ggpmisc)
library(ggpubr)

# Function ----------------------------------------------------------------
GCS_expression_set_count <- function(GCS_df,gtf_df2,bp){
  # count how many GCS falls into each region (upstream, downstream, 5' region, and 3' region)
  # in all genes, genes with high expression, and genes with low expression
  ds =c()
  us =c()
  orf_st =c()
  orf_end =c()
  # all gene features
  for(i in (1:dim(GCS_df)[1])){
    GCS_cor<- GCS_df$pos[i] # genomic position of current GCS
    for (m in (1: dim(gtf_df2)[1])){
      start = gtf_df2$start[m]
      end = gtf_df2$end[m]
      if(start  <= GCS_cor && (start + bp) > GCS_cor && gtf_df2$strand[m] =='+'){
        orf_st <-c(orf_st, GCS_cor)}
      if(start  >= GCS_cor && (start - bp) < GCS_cor && gtf_df2$strand[m] =='-'){
        orf_st <-c(orf_st, GCS_cor)}
      if(end >= GCS_cor && (end - bp) < GCS_cor && gtf_df2$strand[m] =='+'){
        orf_end <-c(orf_end, GCS_cor)}
      if(end <= GCS_cor && (end + bp) > GCS_cor && gtf_df2$strand[m] =='-'){
        orf_end <-c(orf_end, GCS_cor)}
      if(start > GCS_cor && (start - bp) < GCS_cor && gtf_df2$strand[m] =='+'){
        us <-c(us, GCS_cor)}
      if(start < GCS_cor && (start + bp) > GCS_cor && gtf_df2$strand[m] =='-'){
        us <-c(us, GCS_cor)}
      if(end < GCS_cor && (end+bp) > GCS_cor && gtf_df2$strand[m] =='+'){
        ds <-c(ds, GCS_cor)}
      if(end > GCS_cor && (end-bp) < GCS_cor && gtf_df2$strand[m] =='-'){
        ds <-c(ds, GCS_cor)}
    }
  }
  print(paste(('No of GCS in us, 5\', 3\', and ds region'), length(us),length(orf_st),length(orf_end),length(ds)))
  ct_list = list('upstream'=us,'orf_st'= orf_st,'orf_end'=orf_end,'ds'=ds)
  return(ct_list)
}

GCS_expression_strength <- function(GCS_df,bp){
  high_list=GCS_expression_set_count (GCS_df,high_ex,bp)
  low_list=GCS_expression_set_count (GCS_df,low_ex,bp)
  all_list=GCS_expression_set_count (GCS_df,all_ex,bp)

  ##### reshape GCS data frame and remove outliers
  ct = dim(GCS_df)[1]  #number of total GCS   
  GCS_LEVO<-GCS_df %>%
    melt(id.vars=c('pos'),measure.vars =c('strength_rep1','strength_rep2','strength_rep3'),
         variable.name = 'rep',value.name ='strength')

  ##### Combine the count info into a data frame and normalize the strength data
  df_us = GCS_LEVO[which(GCS_LEVO$pos %in% all_list$upstream),] %>% mutate(label = 'upstream',set = 'all_gene')
  df_orf_st = GCS_LEVO[which(GCS_LEVO$pos %in% all_list$orf_st),] %>% mutate(label = '5\'',set = 'all_gene')
  df_orf_end = GCS_LEVO[which(GCS_LEVO$pos %in% all_list$orf_end),]%>%  mutate(label = '3\'',set = 'all_gene')
  df_ds = GCS_LEVO[which(GCS_LEVO$pos %in% all_list$ds),]%>% mutate(label = 'downstream',set = 'all_gene')
  
  df_hius = GCS_LEVO[which(GCS_LEVO$pos %in% high_list$upstream),]%>%
    mutate(label = 'upstream',set = 'high_expression')
  df_hiorf_st = GCS_LEVO[which(GCS_LEVO$pos %in% high_list$orf_st),]%>%
    mutate(label = '5\'',set = 'high_expression')
  df_hiorf_end = GCS_LEVO[which(GCS_LEVO$pos %in%high_list$orf_end),]%>%
    mutate(label = '3\'',set = 'high_expression')
  df_hids = GCS_LEVO[which(GCS_LEVO$pos %in% high_list$ds),]%>%
    mutate(label = 'downstream',set = 'high_expression')
  
  df_lous = GCS_LEVO[which(GCS_LEVO$pos %in% low_list$upstream),]%>%
    mutate(label = 'upstream',set = 'low_expression')
  df_loorf_st = GCS_LEVO[which(GCS_LEVO$pos %in% low_list$orf_st),]%>%
    mutate(label = '5\'',set = 'low_expression')
  df_loorf_end = GCS_LEVO[which(GCS_LEVO$pos %in% low_list$orf_end),]%>%
    mutate(label = '3\'',set = 'low_expression')
  df_lods = GCS_LEVO[which(GCS_LEVO$pos %in% low_list$ds),]%>%
    mutate(label = 'downstream',set = 'low_expression')
  
  df_combined = rbind(df_us,df_orf_st,df_orf_end,df_ds, df_hius,df_hiorf_st,
                      df_hiorf_end,df_hids,df_lous, df_loorf_st, df_loorf_end, df_lods)
  
  info_df = data.frame( count = c(dim(all_ex)[1],dim(high_ex)[1],dim(low_ex)[1]),
                        set = c('all_gene', 'high_expression','low_expression'))
  
  df_combined_strength <-df_combined %>%
    group_by(rep,label,set) %>%
    summarise(total_strength = sum(strength))%>%
    merge(info_df)
  
  df_norm<- GCS_LEVO %>% group_by(rep) %>%
    filter(strength<1e3)%>%
    summarise(strength_sum = sum(strength))%>%
    merge( df_combined_strength)%>%
    mutate(strength_norm = total_strength/(ct*count*300/4642893*(strength_sum/ct))) #ct: No. of cleavge sites, E[No. of cleavage]* Ang stength/site
  return(df_norm)
}


GCS_expression_number <- function(GCS_df,bp){
  high_list=GCS_expression_set_count (GCS_df,high_ex,bp)
  low_list=GCS_expression_set_count (GCS_df,low_ex,bp)
  all_list=GCS_expression_set_count (GCS_df,all_ex,bp)

  # Combine the count info into a dataframe and normalize the count data
  all<-length(all_list$upstream)+length(all_list$orf_st)+length(all_list$orf_end)+length(all_list$ds)
  high_all<-length(high_list$upstream)+length(high_list$orf_st)+length(high_list$orf_end)+length(high_list$ds)
  low_all<-length(low_list$upstream)+length(low_list$orf_st)+length(low_list$orf_end)+length(low_list$ds)
  
  ct = dim(GCS_df)[1]
  
  df<- data.frame(set = rep(c('all_gene', 'high_expression','low_expression'), each =4), 
                  position = (rep(c('upstream','5\'','3\'','downstream'), 3)),
                  count = c(length(all_list$upstream)/(ct*dim(all_ex)[1]*300/4642893), length(all_list$orf_st)/(ct*dim(all_ex)[1]*300/4642893),length(all_list$orf_end)/(ct*dim(all_ex)[1]*300/4642893), length(all_list$ds)/(ct*dim(all_ex)[1]*300/4642893), 
                            length(high_list$upstream)/(ct*dim(high_ex)[1]*300/4642893), length(high_list$orf_st)/(ct*dim(high_ex)[1]*300/4642893), length(high_list$orf_end)/(ct*dim(high_ex)[1]*300/4642893), length(high_list$ds)/(ct*dim(high_ex)[1]*300/4642893),
                            length(low_list$upstream)/(ct*dim(low_ex)[1]*300/4642893), length(low_list$orf_st)/(ct*dim(low_ex)[1]*300/4642893), length(low_list$orf_end)/(ct*dim(low_ex)[1]*300/4642893), length(low_list$ds)/(ct*dim(low_ex)[1]*300/4642893)),
                  count_raw = c(length(all_list$upstream), length(all_list$orf_st), length(all_list$orf_end), length(all_list$ds), 
                                length(high_list$upstream), length(high_list$orf_st), length(high_list$orf_end), length(high_list$ds), 
                                length(low_list$upstream), length(low_list$orf_st), length(low_list$orf_end), length(low_list$ds))
  )
  return(df)
}

# Input:  -------------------------------------------------------------------
## read annotation file
.anoPath <- '/data/REFERENCE/GCF_000005845.2_ASM584v2_genomic_insert_hierexon.gffread.gtf'
.colselect = c("gene_id","gene_name",'width','start','end','strand')
gtf_df<-readAnnotaion(AnnotationPath =.anoPath,col_to_select=.colselect)


##read GCS file
### format: /pos/tMOXI_strength_rep1/tMOXI_strength_rep2/tMOXI_strength_rep3
.gcspath <- '/data/GCS_calling/GCS_extended_version/'
.suffix <- '.txt'
drug = c('MOXI','LEVO','GEMI','CIP','NOR')
FQ_strength_list <- lapply(drug,readGCSStrength, GCSPath = .gcspath, suffix  =.suffix,move_first = FALSE,col_to_select=c(1,24,17,20,23))
names(FQ_strength_list)=drug 


## Read in sample expression data
.RNAseqpath <- 'data/gene_exp/COUNTS/RAW'
.condition = c('treatment') #sub folder name
.drug = c('negative') #sub folder name
neg_counts <- getRNAseqcounts(RNAseqPath=.RNAseqpath,condition=.condition,drug=.drug) # Load in raw count data 
colnames(neg_counts) <- c('gene_id','rep1','rep2','rep3')
neg_counts <- neg_counts %>% left_join(gtf_df, by='gene_id') #concatenate with annotation

### normalized by gene length (transcript/kb per million reads)
neg_30min_TPM <-neg_counts %>% 
  transform (rep1_norm = (rep1+1) /(width/1000),rep2_norm = (rep2+1)/(width/1000), rep3_norm = (rep3+1) /(width/1000))%>% 
  transform(TPM_rep1 =rep1_norm/(sum(rep1_norm)/1000000), TPM_rep2 =rep2_norm/(sum(rep2_norm)/1000000),TPM_rep3 =rep3_norm/(sum(rep3_norm)/1000000) )

# Bin the gene expression level into low and high  
high_ex<-neg_30min_TPM %>% 
  transform(TPM_avg = (TPM_rep1+TPM_rep2+TPM_rep3)/3)%>%
  filter(TPM_avg>400)%>% # 382 genes
  transform(label= 'high')
high_gene_ct = dim(high_ex)[1]
print(paste('No. of highly expressed genes: ', high_gene_ct))

low_ex<-neg_30min_TPM %>% 
  transform(TPM_avg = (TPM_rep1+TPM_rep2+TPM_rep3)/3)%>%
  filter(TPM_avg<1.5)%>% # 404 genes
  transform(label= 'low')
low_gene_ct = dim(low_ex)[1]
print(paste('No. of low expressed genes: ', low_gene_ct))

all_ex<- neg_30min_TPM %>% 
  transform(TPM_avg = (TPM_rep1+TPM_rep2+TPM_rep3)/3)%>% 
  transform(label= 'all')
all_gene_ct = dim(all_ex)[1]



########
#SAVE the number of GCS to excel
e1 = GCS_expression_strength(FQ_strength_list$MOXI,300)
e2 = GCS_expression_strength(FQ_strength_list$LEVO,300)
e3 = GCS_expression_strength( FQ_strength_list$GEMI,300)
e4 = GCS_expression_strength(FQ_strength_list$CIP,300)
e5 = GCS_expression_strength(FQ_strength_list$NOR,300)
write.csv(rbind(e1,e2,e3,e4,e5), "GCS_expression_strength_sets.csv")


df1<-GCS_expression_number(FQ_strength_list$MOXI,300)
df2 <-GCS_expression_number(FQ_strength_list$LEVO,300)
df3<-GCS_expression_number(FQ_strength_list$GEMI,300)
df4 <-GCS_expression_number(FQ_strength_list$CIP,300)
df5 <-GCS_expression_number(FQ_strength_list$NOR,300)
write.csv(rbind(df1,df2,df3,df4,df5), "GCS_expression_number_sets.csv")


# count plot---------------------------------------------------------
ct_plot <- function(df,title){
  s1<- ggplot(df,aes(x=set, y=count,fill = fct_inorder(position)))+
    geom_bar(stat= "identity", position = position_dodge(width = 0.9),color='black')+
    xlab( '')+
    ylab( "")+
    scale_fill_brewer(palette="GrBG")+
    theme_classic()+
    ggtitle(title)+
    theme(plot.title = element_text(size =10), legend.title = element_blank())
  return(s1)
}

c1 <- ct_plot(df1,'MOXI')
c2 <- ct_plot(df2,'LEVO')
c3 <- ct_plot(df3,'GEMI')
c4 <- ct_plot(df4,'CIP')
c5 <- ct_plot(df5,'NOR')

c_sum<-ggarrange(c1,c2,c3,c4,c5, ncol=1, nrow=5, common.legend = TRUE, legend="bottom")
ggsave("GCS_sets_300bp_count.pdf", units="in", c_sum, width=3.5, height=9, dpi=1200)


# strength plot  ---------------------------------------------------------
str_plot <- function(df,title){
  df$label = factor(df$label,levels = c('upstream','5\'','3\'','downstream'))
  s1<-ggbarplot(df,x='set', y='strength_norm', add =  c("mean_se"),
                fill = 'label',alpha=0.8,width =0.9,position = position_dodge(0.9))+
    xlab('')+ylab('')+ggtitle(title)+
    scale_y_continuous(breaks=seq(0,1,by=0.5))+
    scale_fill_brewer(palette="GrBG")+
    theme_classic()+
    theme(plot.title = element_text(size =10),legend.title = element_blank())
  
  return(s1)
}
s1 <- str_plot(e1,'MOXI')
s2 <- str_plot(e2,'LEVO')
s3 <- str_plot(e3,'GEMI')
s4 <- str_plot(e4,'CIP')
s5 <- str_plot(e5,'NOR')

s_sum<-ggarrange(s1,s2,s3,s4,s5, ncol=1, nrow=5, common.legend = TRUE, legend="bottom")
ggsave("GCS_sets_300bp_expression.pdf", units="in", s_sum, width=3.5, height=7, dpi=1200)


# s2<-ggbarplot(e2,x='set', y='strength_norm', add =  c("mean_se"),
#               fill = 'label',alpha=0.8,width =0.9,position = position_dodge(0.9))+
#   xlab('')+ylab('')+ggtitle("LEVO")+
#   scale_y_continuous(breaks=seq(0,1,by=0.5))+
#   scale_fill_brewer(palette="GrBG")+
#   theme_classic()+
#   theme(plot.title = element_text(size =10),legend.title = element_blank())

# stats, strength ---------------------------------------------------------
ano_test <- function(e1){
  df = e1 %>%
    filter(set=='all_gene')
  print('Performing One-way anova,all_gene...')
  Aov1<- aov(total_strength ~ label, data = df)
  print(summary(Aov1))
  plot(Aov1, 1)
  print('Performing Tukey...')
  print(TukeyHSD(Aov1))
  
  df =e1 %>%
    filter(set=='high_expression')
  print('Performing One-way anova,high_expression...')
  Aov1<- aov(total_strength ~ label, data = df)
  print(summary(Aov1))
  print(TukeyHSD(Aov1))
  
  df = e1 %>%
    filter(set=='low_expression')
  print('Performing One-way anova,low_expression...')
  Aov1<- aov(total_strength ~ label, data = df)
  print(summary(Aov1))
  print(TukeyHSD(Aov1))
}

ano_test(e1)
ano_test(e2)
ano_test(e3)
ano_test(e4)
ano_test(e5)





###### binomial test --------------------------------------------------------------
# LEVO GCS, total of 4170 genes
binom.test(m, sum_ct, p = prob, alternative = c("two.sided"))







