library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(GenomicRanges)
library(dplyr)


Flag<-read.table("data/GCS_calling/coverage_depth/LEVO/Flag1.txt")
noFlag<-read.table("data/GCS_calling/coverage_depth/LEVO/Flagless1.txt")
df<-data.frame(Flag = Flag$Flag1,Flagless =noFlag$Flagless1, pos = Flag$Pos) 

mu<-df[3912487:3913087,] %>% #600bp windwos
  melt(id.vars = ('pos'),variable.name = 'Sample', value.name = 'coverage')

nuoN<-df[2389880:2390480,]%>%
  melt(id.vars = ('pos'),variable.name = 'Sample', value.name = 'coverage')
mobA<-df[4042564:4043164,]%>%
  melt(id.vars = ('pos'),variable.name = 'Sample', value.name = 'coverage')
trkH<-df[4035475:4036075,]%>%
  melt(id.vars = ('pos'),variable.name = 'Sample', value.name = 'coverage')


# plot --------------------------------------------------------------------
s1 <- ggplot(mu)+
  geom_line(aes(x=pos, y = coverage,group =Sample,color = Sample))+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  xlab('[3,912,487, 3,913,087]')+
  ylab('Coverage')+
  ggtitle('Mu')+
  geom_rect(aes(xmin = 3912641, xmax=3912766, ymin= -1800, ymax=-300), fill = 'green',alpha = 0.4)+
  theme(axis.title.x=element_blank(),legend.position = 'none',legend.title = element_blank(),plot.title = element_text(vjust = - 5,hjust =0.1))
  s1
  
s2<-ggplot(mobA)+ #4,041,415 <- 4,041,999
  geom_line(aes(x=pos, y = coverage,group =Sample,color = Sample))+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  xlab('')+
  ylab('Coverage')+
  ggtitle('mobA')+ #4042864
  geom_rect(aes(xmin = 4042962, xmax=4043036, ymin= 0, ymax=150), fill = 'green',alpha = 0.4)+
  theme(axis.title.x=element_blank(),legend.position = 'none',legend.title = element_blank(),plot.title = element_text(vjust = - 5,hjust =0.1,face="italic"))

s3<-ggplot(nuoN)+ #2,390,048<- 2,391,505
  geom_line(aes(x=pos, y = coverage,group =Sample,color = Sample))+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  xlab('')+
  ylim(0,6000)+
  ggtitle('nuoN')+
  geom_rect(aes(xmin = 2390126, xmax=2390234, ymin= 0, ymax=120), fill = 'green',alpha = 0.4)+
  ylab('Coverage')+
  theme(axis.title.x=element_blank(),legend.position = 'none',legend.title = element_blank(),plot.title = element_text(vjust = - 5,hjust =0.1,face="italic"))


s4<-ggplot(trkH)+ 
  geom_line(aes(x=pos, y = coverage,group =Sample,color = Sample))+
  scale_color_brewer(palette = "Set1")+
  theme_classic()+
  xlab('')+
  ylim(0,6000)+
  ggtitle('trkH')+
  ylab('Coverage')+
  geom_rect(aes(xmin = 4035723, xmax=4035828, ymin= 0, ymax=120), fill = 'green',alpha = 0.4)+
  theme(axis.title.x=element_blank(),legend.position = 'none',legend.title = element_blank(),plot.title = element_text(vjust = - 5,hjust =0.1,face="italic"))

ss_sum <- grid.arrange(s1,s2,s3,s4, nrow =4)  
ggsave("coverage.pdf", units="in", s_sum, width=3.5, height=6, dpi=1200)

######################

# Coverage depth plot (Fig1) ----------------------------------------------------
#Mu 
Flag<-c(58828,68739,58798,58784,33513,33494,33612,34178,66980,69019,69076,66571)#4 b4, 4after
noFlag<-c(1105,1106,1101,1088,850,861,868,865,1058,1119,1153,1155)

m1<-barplot(Flag,ylim=c(0,69100))
m2<-barplot(noFlag,ylim=c(0,69100))
#barplot(noFlag,ylim=c(0,1500))


#mobA

Flag<-c(5431,6121,6667,6739,5165,4734,4715,4723,5897,6375,6422,6446)#4 b4, 4after
noFlag<-c(1097,1088,1104,1096,1041,1038,1037,1045,1038,1028,1033,1047)
m3<-barplot(Flag,ylim=c(0,6500))
m4<-barplot(noFlag,ylim=c(0,6500))
#barplot(noFlag,ylim=c(0,1500))

s_sum <- grid.arrange(s1,s2,s3,s4,s_sum <- grid.arrange(s1,s2,s3,s4, nrow =4)  , nrow =2)  

######################
#Strain info
library(Biostrings)  
muori = readDNAStringSet('/Volumes/TJC/GCS_calling_files/Mu_ori_mu_insert_MG1655.fa')
Biostrings::letterFrequency(muori, letters ='ATCG',as.prob = TRUE)
alphabetFrequency(muori,baseOnly=TRUE, as.prob=TRUE)
