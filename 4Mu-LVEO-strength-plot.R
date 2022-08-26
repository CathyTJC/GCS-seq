## Mu strength and LEVO-Mu comparison
library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(GenomicRanges)
library(ggpmisc)
library(ggpubr)

# 4mu strain Mu strength --------------------------------------------------
# extract Mu site strength (1567679, No.10, yddy and yddW; 1806970, ds bfkB;3913663 glmS; 3925233 No. 54 ds rsmG)
setwd("/Volumes/TJC/GCS_calling_files/")

Mu4_LEVO_raw <- read.table("HF_GCS_strength_0408/4-Mu/4Mu-LEVO.txt",sep = ',',header = TRUE)
Mu4_LEVO <- read.table("HF_GCS_strength_0408/4-Mu/4Mu-LEVO.txt",sep = ',',header = TRUE) %>%
  filter(pos %in% c(1567679,1806970,3913663,3925233))%>%
  select(pos, rep1=strength_rep1, rep2=strength_rep2, rep3=strength_rep3)%>%
  gather(key='rep', value = 'Strength', rep1:rep3,-pos) %>%
  arrange(pos) %>%
  transform(Strain = 'Mu_ori_ter', Position = as.factor(rep(c('1567679','1806970','3913663','3925233'), each =3)))


## Mu expression varies at different location
Mu4_LEVO_strength_plot <-ggplot(Mu4_LEVO,aes(x=Position, y = Strength,fill=Position))+
  geom_boxplot()+
  scale_fill_brewer(palette="RdBu") + 
  theme_classic()+
  xlab('Mu GCS Position')+
  ylab('Cleavage Strength')
print(Mu4_LEVO_strength_plot)

ggsave('Mu4-LEVO.png',Mu4_LEVO_strength_plot,units="in", width=5, height=3, dpi=1200)



# LEVO-4LEVO comparison ---------------------------------------------------
setwd("/Volumes/TJC/GCS_calling_files/LEVO_compare_strength/")

LEVO <-read.table('LEVO.txt', sep = ',',header = TRUE)%>%
  select('pos','strength_rep1','strength_rep2','strength_rep3')%>%
  replace(., is.na(.), 0) %>%
  filter(strength_rep1<1e9,strength_rep2<1e9,strength_rep3<1e9)

LEVO4Mu <- read.table("4Mu-LEVO-corrected.txt",sep = ',', header = TRUE)%>%
  select('mod_GCS','strength_rep1','strength_rep2','strength_rep3')%>%
  replace(., is.na(.), 0)%>%
  filter(strength_rep1<1e9,strength_rep2<1e9,strength_rep3<1e9)
colnames(LEVO4Mu) = c('pos','strength_rep1','strength_rep2','strength_rep3')

s1 = LEVO %>%
  filter(strength_rep1<1e9,strength_rep2<1e9,strength_rep3<1e9)%>%
  ggplot(aes(strength_rep1,strength_rep2))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO-rep1')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('LEVO-rep2')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

s2 = LEVO %>%
  filter(strength_rep1<1e9,strength_rep2<1e9,strength_rep3<1e9)%>%
  ggplot(aes(strength_rep1,strength_rep3))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO-rep1')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('LEVO-rep3')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

s3 = LEVO %>%
  filter(strength_rep1<1e9,strength_rep2<1e9,strength_rep3<1e9)%>%
  ggplot(aes(strength_rep2,strength_rep3))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO-rep2')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('LEVO-rep3')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()


p1 = LEVO4Mu %>%
  filter(strength_rep1<1e9,strength_rep2<1e9,strength_rep3<1e9)%>%
  ggplot(aes(strength_rep1,strength_rep2))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('4Mu-LEVO-rep1')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep2')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

p2 = LEVO4Mu %>%
  filter(strength_rep1<1e9,strength_rep2<1e9,strength_rep3<1e9)%>%
  ggplot(aes(strength_rep1,strength_rep3))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('4Mu-LEVO-rep1')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep3')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

p3 = LEVO4Mu %>%
  filter(strength_rep1<1e9,strength_rep2<1e9,strength_rep3<1e9)%>%
  ggplot(aes(strength_rep2,strength_rep3))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('4Mu-LEVO-rep2')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep3')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

df_merge1 = LEVO %>% select('pos','strength_rep1') %>%
  merge(LEVO4Mu,by = 'pos')

df_merge2 = LEVO %>% select('pos','strength_rep2') %>%
  merge(LEVO4Mu,by = 'pos')

df_merge3 = LEVO %>% select('pos','strength_rep3') %>%
  merge(LEVO4Mu,by = 'pos')


g1 = df_merge1 %>%
  ggplot(aes(strength_rep1.x,strength_rep1.y))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep1')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep1')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

g2 = df_merge1 %>%
  ggplot(aes(strength_rep1.x,strength_rep2))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep1')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep2')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()


g3 = df_merge1 %>%
  ggplot(aes(strength_rep1.x,strength_rep3))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep1')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep3')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

g4 = df_merge2 %>%
  ggplot(aes(strength_rep2.x,strength_rep1))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep2')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep1')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

g5 = df_merge2 %>%
  ggplot(aes(strength_rep2.x,strength_rep2.y))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep2')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep2')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

g6 = df_merge2 %>%
  ggplot(aes(strength_rep2.x,strength_rep3))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep2')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep3')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()


g7 = df_merge3 %>%
  ggplot(aes(strength_rep3.x,strength_rep1))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep3')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep1')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

g8 = df_merge3 %>%
  ggplot(aes(strength_rep3.x,strength_rep2))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep3')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep2')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

g9 = df_merge3 %>%
  ggplot(aes(strength_rep3.x,strength_rep3.y))+geom_point()+
  geom_abline(colour = "brown")+
  xlab('LEVO_rep3')+
  xlim(0,8)+
  ylim(0,8)+
  ylab('4Mu-LEVO-rep3')+
  stat_poly_eq(formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(label = paste(..rr.label..)), 
               parse = TRUE)+
  theme_classic()

s_sum<-ggarrange(s1,s2,s3,p1,p2,p3,g1,g2,g3,g4,g5,g6,g7,g8,g9, ncol=3, nrow=5)
ggsave("Mu-Mu-compare.png", units="in", s_sum, width=7, height=8, dpi=1200)


# ggsave(avg_strength_plot, '4muvslevo strength', units="in",  width=3.5, height=3.5, dpi=1200)




  
##################################################
# count<-c(7,10,0,8,2)
# Drug <-c('LEVO','MOXI','NOR','GEMI','CIP')
# count_50<-c(38,48,0,34,15)
# df_ct<-data.frame(count, Drug)
# 
# s1<-ggplot(df_ct,aes(x=Drug, y=count,fill=Drug))+
#   geom_bar(stat = "identity") +
#   scale_fill_brewer(palette="Blues")+
#   theme_classic()
# 
# 
# s2<-ggplot(df_ct,aes(x=Drug, y=count_50,fill=Drug))+
#   geom_bar(stat = "identity") +
#   scale_fill_brewer(palette="Blues")+
#   theme_classic()
# 
# s_sum <- grid.arrange(s1,s2, nrow =1)  
# ggsave("GCS_ct.png", units="in", s_sum, width=10, height=8, dpi=600)

