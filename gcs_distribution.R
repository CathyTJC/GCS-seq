############
# This script plot the GCS distribution across the chromosome and perform statistical tests
# @Author Cathy Tang
############

###### --------------------------------------------------------------
# Library import
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(ggforce)
library(car)
library(userfriendlyscience)
library(ggpubr)

gcsCount<-function(gcs,bin_info){
  gcs_ct = gcs[,1]
  bin_ct<-bin_info %>%
    rowwise()%>%
    mutate(count = sum(gcs_ct>=start & gcs_ct<=end))
  return(bin_ct)
}

gcsCount_strength<-function(gcs,bin_info){
  gcs_ct = gcs$pos
  gcs_strength1 = gcs$strength_rep1
  gcs_strength2 = gcs$strength_rep2
  gcs_strength3 = gcs$strength_rep3
  gcs_strength = gcs$avg_strength
  bin_ct<-bin_info %>%
    rowwise()%>%
    mutate(count1 = sum(gcs_strength1[which(gcs_ct>=start &  gcs_ct<=end &  gcs_strength<1e10)]),
           count2 = sum(gcs_strength2[which(gcs_ct>=start &  gcs_ct<=end &  gcs_strength<1e10)]),
           count3 = sum(gcs_strength3[which(gcs_ct>=start &  gcs_ct<=end &  gcs_strength<1e10)]))%>%
    melt(id.vars=c('bin','start','end'),measure.vars =c('count1','count2','count3'),
         variable.name = 'rep',value.name ='count')
  
  return(bin_ct)
}

#normalize by MD length (No. of distinct GCS)
MD_normalize<-function(gcs){
  gcs$length = gcs$end-gcs$start+1
  gcs<- gcs %>%
    select(c('bin', 'count','length')) %>%
    group_by (bin) %>%
    summarise(count = sum(count), length = sum(length))
  gcs$prob = gcs$length/4642893
  gcs$norm_count = gcs$count/(sum(gcs$count)*gcs$length/4642893)
  return (gcs)
}

MD_strength_normalize<-function(gcs,ct_df){
  total_count = sum(ct_df$count)
  gcs$length = gcs$end-gcs$start+1
  gcs<- gcs %>%
    select(c('bin', 'count','length','rep')) %>%
    group_by (bin,rep) %>%
    summarise(count = sum(count), length = sum(length))
  gcs$prob = gcs$length/4642893
  
  tot_strength <-gcs %>%
    group_by (rep)%>%
    summarise(total_strength = sum(count))
  gcs<-merge(gcs,tot_strength)
  
  gcs$norm_strength = gcs$count/(total_count*(gcs$length/4642893)*gcs$total_strength/total_count)
  return (gcs)
}

binom<-function(df,p){
  ct<-df$count
  sum_ct<-sum(ct)
  print(sum_ct)
  binnum = 0
  for (m in 1: length(ct)){
    count = ct[m]
    pval = p[m]
    binnum = binnum+1
    print(paste0('Stats for Bin ', binnum))
    stas<-binom.test(count , sum_ct, pval, alternative = c("two.sided"))
    print(stas)
  }
}

###### input --------------------------------------------------------------
# Read in GCS files:LEVO, GEMI, LEVO, MOXI, NOR HF_GCSs
dataDir ='/data/GCS_calling/GCS_extended_version/'
dfFiles <- dir(file.path(getwd(), dataDir), pattern="*.txt", full.name=T) 
names(dfFiles) <- gsub(".txt","",basename(dfFiles))

CIP<-read.table(dfFiles[1],sep = ',',header = TRUE)[,c(1,24,17,20,23)]
GEMI<-read.table(dfFiles[2],sep = ',',header = TRUE)[,c(1,24,17,20,23)]
LEVO<-read.table(dfFiles[3],sep = ',',header = TRUE)[,c(1,24,17,20,23)]
MOXI<-read.table(dfFiles[4],sep = ',',header = TRUE)[,c(1,24,17,20,23)]
NOR<-read.table(dfFiles[5],sep = ',',header = TRUE)[,c(1,24,17,20,23)]


# Bin information construction: A:J 
bin = 464289 
start<-seq(1,464290+8*bin,bin)
end<-seq(464289,4642893,bin)
end[10]<-4642893
Bin<-LETTERS[1:10]
bin_info<-data.frame(bin = factor(Bin), start=start, end=end)

# MD information construction
Bin<-c('Ori','NS-R','Right','Ter','Left','NS-L','Ori')
level <-c('NS-R','Right','Ter','Left','NS-L','Ori')
start <-c(1,46418,603416,1206831,2181577,2877825,3759739)
end<-c(46417,603415,1206830,2181576,2877824,3759738,4642893)
MD_info<-data.frame(bin = factor(Bin,levels = level), start=start, end=end)

######  number of distinct GCSs,10bin--------------------------------------------------------------
# plot GCS distribution (number of distinct GCSs)
LEVO_ct <-gcsCount(LEVO,bin_info)
CIP_ct <-gcsCount(CIP,bin_info)
GEMI_ct <-gcsCount(GEMI,bin_info)
MOXI_ct <-gcsCount(MOXI,bin_info)
NOR_ct <-gcsCount(NOR,bin_info)

LEVO_pt<-ggplot(LEVO_ct,aes(x=bin, y=count))+
  geom_bar(fill="skyblue2",color ='Black',alpha =0.8,stat = "identity")+
  theme_bw()+
  xlab('')+
  ylim(0,1500)+
  ylab('')+
  ggtitle("LEVO")+
  theme(plot.title = element_text(size =10))

CIP_pt<-ggplot(CIP_ct,aes(x=bin, y=count))+
  geom_bar(fill="skyblue2",color ='Black',alpha =0.8,stat = "identity")+
  theme_bw()+
  xlab('')+
  ylim(0,1500)+
  ylab('')+
  ggtitle("CIP")+
  theme(plot.title = element_text(size =10))

GEMI_pt<-ggplot(GEMI_ct,aes(x=bin, y=count))+
  geom_bar(fill="skyblue2",color ='Black',alpha =0.8,stat = "identity")+
  theme_bw()+
  xlab('')+
  ylab('Number of distinct GCSs')+
  ylim(0,1500)+
  ggtitle("GEMI")+
  theme(plot.title = element_text(size =10))

MOXI_pt<-ggplot(MOXI_ct,aes(x=bin, y=count))+
  geom_bar(fill="skyblue2",color ='Black',alpha =0.8,stat = "identity")+
  theme_bw()+
  xlab('')+
  ylab('')+
  ggtitle("MOXI")+
  theme(plot.title = element_text(size =10))

NOR_pt<-ggplot(NOR_ct,aes(x=bin, y=count))+
  geom_bar(fill="skyblue2",color ='Black',alpha =0.8,stat = "identity")+
  theme_bw()+
  ylim(0,1500)+
  xlab('')+
  ylab('')+
  ggtitle("NOR")+
  theme(plot.title = element_text(size =10))


NOR_pt_zm<-ggplot(NOR_ct,aes(x=bin, y=count))+
  geom_bar(fill="skyblue2",color ='Black',alpha =0.8,stat = "identity")+
  theme_bw()+
  xlab('')+
  ylab('')+
  ggtitle("NOR")+
  theme(plot.title = element_text(size =10))


s_sum <- grid.arrange(MOXI_pt,LEVO_pt,GEMI_pt, CIP_pt, NOR_pt, ncol=1) 
ggsave("GCS_distribution_all.png", units="in", s_sum, width=3.5, height=6, dpi=600) #save the image
ggsave("NOR_distribution_zoomin.png", units="in", NOR_pt_zm, width=3.5, height=1.2, dpi=600)


write.csv(rbind(MOXI_ct,LEVO_ct,GEMI_ct,CIP_ct,NOR_ct), "GCSs_distribtuion.csv")

######  strength, 10bin--------------------------------------------------------------
LEVO_strength <-gcsCount_strength(LEVO,bin_info)
CIP_strength <-gcsCount_strength(CIP,bin_info)
GEMI_strength <-gcsCount_strength(GEMI,bin_info)
MOXI_strength <-gcsCount_strength(MOXI,bin_info)
NOR_strength <-gcsCount_strength(NOR,bin_info)


LEVO_pt<-ggbarplot(LEVO_strength,x='bin', y='count', add =  c("mean_se"),
          fill = 'skyblue2',alpha=0.8,width =0.9)+
  xlab('')+
  ylim(0,557)+
  ylab('')+
  ggtitle("LEVO")+
  theme_bw()+
  theme(plot.title = element_text(size =10))

MOXI_pt<-ggbarplot(MOXI_strength,x='bin', y='count', add =  c("mean_se"),
                   fill = 'skyblue2',alpha=0.8,width =0.9)+
  xlab('')+
  ylim(0,557)+
  ylab('')+
  ggtitle("MOXI")+
  theme_bw()+
  theme(plot.title = element_text(size =10))

CIP_pt<-ggbarplot(CIP_strength,x='bin', y='count', add =  c("mean_se"),
                   fill = 'skyblue2',alpha=0.8,width =0.9)+
  xlab('')+
  ylim(0,557)+
  ylab('')+
  ggtitle("CIP")+
  theme_bw()+
  theme(plot.title = element_text(size =10))
  

GEMI_pt<-ggbarplot(GEMI_strength,x='bin', y='count', add =  c("mean_se"),
                   fill = 'skyblue2',alpha=0.8,width =0.9)+
  xlab('')+
  ylim(0,557)+
  ylab('')+
  ggtitle("GEMI")+
  ylab('Cleavage Strength')+
  theme_bw()+
  theme(plot.title = element_text(size =10))

NOR_pt<-ggbarplot(NOR_strength,x='bin', y='count', add =  c("mean_se"),
                   fill = 'skyblue2',alpha=0.8,width=0.9)+
  xlab('')+
  ylim(0,557)+
  ylab('')+
  ggtitle("NOR")+
  theme_bw()+
  theme(plot.title = element_text(size =10))

NOR_pt_zm<-ggbarplot(NOR_strength,x='bin', y='count', add =  c("mean_se"),
                   fill = 'skyblue2',alpha=0.8,width =0.9)+
  xlab('')+
  ylab('')+
  theme_bw()+
  theme(plot.title = element_text(size =10))
  
  
s_sum <- grid.arrange(MOXI_pt,LEVO_pt,GEMI_pt, CIP_pt, NOR_pt, ncol=1) 
ggsave("GCS_strength_distribution.png", units="in", s_sum, width=3.5, height=6, dpi=600) #save the image
ggsave("GCS_NOR-strength_distribution.png", units="in", NOR_pt_zm, width=3.5, height=1.2, dpi=600) #save the image

write.csv(rbind(MOXI_strength,LEVO_strength,GEMI_strength,CIP_strength,NOR_strength), "GCSs_strength_distribtuion.csv")

# ANOVA, strength test ----------------------------------------------------
ano_test<-function(df){
  # test equal variance
  boxplot(count ~ bin,  data=df)
  # print('levene_Test...')
  # print(leveneTest(count ~ bin, data = df))
  
  ##########
  # Anova & kruskal
  ##########
  print('Performing One-way anova...')
  Aov1<- aov(count ~ bin, data = df)
  print(summary(Aov1))
  plot(Aov1, 1)
  print('Performing Tukey...')
  print(TukeyHSD(Aov1))
  
  print('Performing kruskal.test...')
  print(kruskal.test(count ~ bin, data = df))
        # check normality
        res<-Aov1$residuals
        hist(res, main="Histogram of standardised
  residuals",xlab="Standardised residuals")
        qqnorm(res)
        qqline(res)
}

ano_test(LEVO_strength)
ano_test(MOXI_strength)
ano_test(GEMI_strength)
ano_test(CIP_strength)
ano_test(NOR_strength)






# GCS count, MD as bin ----------------------------------------------------
LEVO_ct <-gcsCount(LEVO,MD_info)
CIP_ct <-gcsCount(CIP,MD_info)
GEMI_ct <-gcsCount(GEMI,MD_info)
MOXI_ct <-gcsCount(MOXI,MD_info)
NOR_ct <-gcsCount(NOR,MD_info)

LEVO_ctn<-MD_normalize(LEVO_ct)
CIP_ctn<-MD_normalize(CIP_ct)
GEMI_ctn <-MD_normalize(GEMI_ct)
MOXI_ctn <-MD_normalize(MOXI_ct)
NOR_ctn <-MD_normalize(NOR_ct)

# ggplot(LEVO_ctn,aes(x=bin, y=norm_count))+
#   geom_bar(fill="skyblue2",color ='Black',alpha =0.8,stat = "identity")+
#   theme_bw()+
#   xlab('')+
#   ylab('')+
#   ggtitle("LEVO")+
#   scale_y_continuous(breaks=seq(0.6,1.5,by=0.3))+
#   theme(plot.title = element_text(size =10))

LEVO_pt<-ggplot(LEVO_ctn, aes(x=bin, y=norm_count)) +
  geom_segment(aes(x=bin, xend=bin, y=1, yend=norm_count),alpha = 0.8,color="blue",size=13) +
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dotted')+
  xlab('')+
  ylab('')+
  ggtitle("LEVO")+
  scale_y_continuous(limits = c(0.9,1.1))+
  theme(plot.title = element_text(size =10))

CIP_pt<-ggplot(CIP_ctn, aes(x=bin, y=norm_count)) +
  geom_segment(aes(x=bin, xend=bin, y=1, yend=norm_count), alpha = 0.8,color="blue",size=13) +
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dotted')+
  xlab('')+
  ylab('')+
  ggtitle("CIP")+
  scale_y_continuous(limits = c(0.9,1.1))+
  theme(plot.title = element_text(size =10))

GEMI_pt<-ggplot(GEMI_ctn, aes(x=bin, y=norm_count)) +
  geom_segment(aes(x=bin, xend=bin, y=1, yend=norm_count),alpha = 0.8,color="blue",size=13) +
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dotted')+
  xlab('')+
  ylab('')+
  ggtitle("GEMI")+
  scale_y_continuous(limits = c(0.9,1.1))+
  theme(plot.title = element_text(size =10))

MOXI_pt<-ggplot(MOXI_ctn, aes(x=bin, y=norm_count)) +
  geom_segment(aes(x=bin, xend=bin, y=1, yend=norm_count), alpha = 0.8,color="blue",size=13) +
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dotted')+
  xlab('')+
  ylab('')+
  ggtitle("MOXI")+
  scale_y_continuous(limits = c(0.9,1.1))+
  theme(plot.title = element_text(size =10))

NOR_pt<-ggplot(NOR_ctn, aes(x=bin, y=norm_count)) +
  geom_segment(aes(x=bin, xend=bin, y=1, yend=norm_count), alpha = 0.8,color="blue",size=13) +
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dotted')+
  xlab('')+
  ylab('')+
  ggtitle("NOR")+
  scale_y_continuous(limits = c(0.9,1.1))+
  theme(plot.title = element_text(size =10))

s_sum2 <- grid.arrange(MOXI_pt,LEVO_pt,GEMI_pt, CIP_pt, NOR_pt, ncol=1) 
ggsave("GCS_distributionMD.png", units="in", s_sum2, width=3.5, height=6, dpi=600) #save the image

write.csv(rbind(MOXI_ctn,LEVO_ctn,GEMI_ctn,CIP_ctn,NOR_ctn), "GCSs_ct_distribtuion_MD.csv")


# statistical analysis ----------------------------------------------------
p1 = rep(0.1,10)
binom(LEVO_ct,p1) 
#Bin 4 p-value = 8.046e-07
#Bin 6 p-value = 0.002501
binom(GEMI_ct,p1) 
#6 p-value = 9.276e-05*
binom(MOXI_ct,p1) 
#4 p-value = 3.593e-08
binom(CIP_ct,p1)
#4:p-value = 0.001056
#6: p-value = 0.002998
# 7: 0.004167
binom(NOR_ct,p1) 
#4 p-value = 0.0006978
#6 :0.000813
#7; 0.001727
#8 p-value = 9.581e-05

#statistical analysis, MD
p_MD =LEVO_ctn$prob
binom(LEVO_ctn,p_MD)
binom(GEMI_ctn,p_MD) 
binom(MOXI_ctn,p_MD) 
binom(CIP_ctn,p_MD)
binom(NOR_ctn,p_MD) 





# Strength plot,  MD -----------------------------------------------------------
LEVO_strength <-gcsCount_strength(LEVO,MD_info)
CIP_strength <-gcsCount_strength(CIP,MD_info)
GEMI_strength <-gcsCount_strength(GEMI,MD_info)
MOXI_strength <-gcsCount_strength(MOXI,MD_info)
NOR_strength <-gcsCount_strength(NOR,MD_info)

LEVO_norm<-MD_strength_normalize(LEVO_strength,LEVO_ctn)
CIP_norm <-MD_strength_normalize(CIP_strength,CIP_ctn)
GEMI_norm <-MD_strength_normalize(GEMI_strength,GEMI_ctn)
MOXI_norm <-MD_strength_normalize(MOXI_strength,MOXI_ctn)
NOR_norm <-MD_strength_normalize(NOR_strength,NOR_ctn)


LEVO_st<- ggplot(LEVO_norm, aes(x=bin, y=norm_strength)) +
  geom_boxplot(color = 'skyblue3')+
  geom_jitter(width = 0.2)+
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dashed')+
  xlab('')+
  ylab('')+
  ggtitle("LEVO")+
  scale_y_continuous(limits = c(0.83,1.18))+
  theme(plot.title = element_text(size =9))

MOXI_st<- ggplot(MOXI_norm, aes(x=bin, y=norm_strength)) +
  geom_boxplot(color = 'skyblue3')+
  geom_jitter(width = 0.2)+
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dashed')+
  xlab('')+
  ylab('')+
  ggtitle("MOXI")+
  scale_y_continuous(limits = c(0.83,1.18))+
  theme(plot.title = element_text(size =9))

CIP_st<- ggplot(CIP_norm, aes(x=bin, y=norm_strength)) +
  geom_boxplot(color = 'skyblue3')+
  geom_jitter(width = 0.2)+
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dashed')+
  xlab('')+
  ylab('')+
  ggtitle("CIP")+
  scale_y_continuous(limits = c(0.83,1.18))+
  theme(plot.title = element_text(size =9))

GEMI_st<- ggplot(GEMI_norm, aes(x=bin, y=norm_strength)) +
  geom_boxplot(color = 'skyblue3')+
  geom_jitter(width = 0.2)+
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dashed')+
  xlab('')+
  ylab('')+
  ggtitle("GEMI")+
  scale_y_continuous(limits = c(0.83,1.18))+
  theme(plot.title = element_text(size =9))

NOR_st<- ggplot(NOR_norm, aes(x=bin, y=norm_strength)) +
  geom_boxplot(color = 'skyblue3')+
  geom_jitter(width = 0.2)+
  theme_classic()+
  geom_hline(yintercept = 1,linetype='dashed')+
  xlab('')+
  ylab('')+
  ggtitle("NOR")+
  scale_y_continuous(limits = c(0.83,1.18))+
  theme(plot.title = element_text(size =9))

# CIP_st<-ggbarplot(CIP_norm,x='bin', y='norm_strength', add =  c("mean_se"),
#                   fill = 'skyblue2',alpha=0.8,width =0.9)+
#   xlab('')+ylab('')+ggtitle("CIP")+
#   scale_y_continuous(breaks=seq(0,1,by=0.5))+
#   theme_bw()+theme(plot.title = element_text(size =10))


s_sum <- grid.arrange(MOXI_st,LEVO_st,GEMI_st, CIP_st, NOR_st, ncol=1) 
ggsave("GCS_strength_MD_all_V2.png", units="in", s_sum, width=3.5, height=8, dpi=600) #save the image

write.csv(rbind(MOXI_norm,LEVO_norm,GEMI_norm,CIP_norm,NOR_norm), "GCSs_strength_distribtuion_MD.csv")

# stats, for Strength in MD -------------------------------------------------------------------

Bin<-c('Ori','NS-R','Right','Ter','Left','NS-L')
onesample_ttest <- function(df,Bin){
    for (reg in Bin){
      print(reg)
      sdf <- df %>% filter(bin==reg)
      stats = t.test(sdf$norm_strength, mu=1)
      print(stats)
      if (stats$p.value < 0.05){
        padj = p.adjust(stats$p.value, method = 'BH', 6)
        print(paste0('padj: ',padj))
      }
    }
  return(padj)
}


onesample_ttest(MOXI_norm,Bin) #Ter t = -15.665, df = 2, p-value = 0.00405 padj0.0243015536692717
onesample_ttest(LEVO_norm,Bin) #Ter t = -5.7298, df = 2, p-value = 0.02913 padj: 0.1748
onesample_ttest(GEMI_norm,Bin) # NSL t = -8.0506, df = 2, p-value = 0.01508 padj: 0.0904
# t = -24.238, df = 2, p-value = 0.001698 right  padj: 0.01018
onesample_ttest(CIP_norm,Bin) #ori t = -5.7706, df = 2, p-value = 0.02874  padj: 0.172

onesample_ttest(NOR_norm,Bin) #Ter t = -7.4308, df = 2, p-value = 0.01763  padj: 0.105
