library(ggplot2)
library(tidyverse)
library(RColorBrewer)
library(reshape2)
library(gridExtra)
library(GenomicRanges)
library(ggpmisc)
library(ggpubr)


# GCS number vs persistence level----------------------------------------------------
Drug <-factor(c('MOXI','LEVO','GEMI','CIP','NOR'),levels = c('MOXI','LEVO','GEMI','CIP','NOR'))
Survival_Fraction<-c(7.63E-04,1.96E-02,3.98E-02,8.53E-02,8.11E-01)
SE<- c(3.04E-04,2.44E-03,1.69E-02,3.10E-02,2.39E-01)
GCS_num <-c(14504,12283,11644,6634,1551)
sum_list = c(4474, 2906, 2520, 1697, 323)
sum_list_956 = c(1093, 830, 582, 670, 246)
df1<-data.frame (Treatment=Drug, GCS_num=GCS_num, Survival_Fraction=Survival_Fraction,SE=SE,strength_sum=sum_list)

sur<-c(1.05E-03,1.68E-02,1.79E-02,2.35E-02,1.05E+00,1.08E-03,1.75E-02,
       2.23E-02,1.13E-01,3.33E-01,1.55E-04,2.44E-02,2.89E-02,1.20E-01,1.05E+00)
df<-data.frame(GCS_num=rep(GCS_num,3), GCS_strength =rep(sum_list,3), Treatment=factor(rep(Drug,3),levels = Drug), survival=sur,N = rep(c('rep1','rep2','rep3'),each=5))
df$log_survival <-log10(df$survival)


#plot the fitting curve
lm1 <-ggplot(df1,mapping =aes(y=Survival_Fraction, x=GCS_num))+
  scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)+
  geom_point(aes(color=Treatment),size = 3)+
  theme_bw()+
  xlim(0,15070)+
  geom_errorbar(aes(x=GCS_num ,ymin=Survival_Fraction-SE, ymax=Survival_Fraction+SE,color=Treatment),size=0.5,width=500)+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  theme(legend.position ='none')+
  stat_smooth(df1,mapping = aes(y=Survival_Fraction, x=GCS_num), method = "lm", col = "black",formula = y~x)+
  stat_poly_eq(data = df1, formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(y=Survival_Fraction, x=GCS_num, label = paste(..rr.label..)), 
               parse = TRUE)+
  xlab('Number of distinct GCSs')+
  ylab ('Survival Fraction')+
  annotation_logticks(sides = "l")

lm2 <-ggplot(df1,mapping =aes(y=Survival_Fraction, x=strength_sum))+
  scale_color_brewer(type = 'div', palette = 'Set1', direction = 1)+
  geom_point(aes(color=Treatment),size = 3)+
  theme_bw()+
  geom_errorbar(aes(x=strength_sum ,ymin=Survival_Fraction-SE, ymax=Survival_Fraction+SE,color=Treatment),size=0.5,width=50)+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  stat_smooth(df1, mapping = aes(y=Survival_Fraction, x=strength_sum), method = "lm", col = "black",formula = y~x)+
  stat_poly_eq(data = df1, formula = y~x, label.y = 0.95, label.x = 0.98,
               aes(y=Survival_Fraction, x=strength_sum, label = paste(..rr.label..)), 
               parse = TRUE)+
  xlab('Cleavage Strength')+
  ylab ('Survival Fraction')+
  annotation_logticks(sides = "l")


lm_sum = ggarrange(lm1,lm2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
ggsave("Persistence_GCS_num.pdf", units="in", lm_sum, width=7, height=3.5, dpi=1200)


fit1 <- lm(log(Survival_Fraction)~ strength_sum, data=df1)
summary(fit1) 
AIC(fit1) 
BIC(fit1) 



# bar plot of nubmer of GCS ----------------------------------------------------
p2<-ggplot(df1, aes(x=Treatment, y=GCS_num)) +
  geom_segment(aes(x=Treatment, xend=Treatment, y=0, yend=GCS_num), color="skyblue") +
  geom_point(color="blue", size=4, alpha=0.6) +
  theme_light() +
  coord_flip() +
  theme_classic()+
  xlab('')+
  theme(text = element_text(size = 14)) +  
  expand_limits(y = c(0, 15000))+
  ylab('Number of distinct GCSs')

ggsave("GCS_persistance_bar.pdf",p2,width=4, height=3, units = "in", dpi = 600)




