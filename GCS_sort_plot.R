################
# Sort GCS based on the cleavage strength and save the file of sites and strength for motif construction
# Further, retain the top sites which are identified as real GCSs for motif construction 
# Also heatmap and barplots showng the  top strengths sites for each FQ treatment 
# next step in pipeline for motif construction: pssm_construction.py
###### 
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)
library(tidyr)
library(data.table)
library(tidyverse)

save_df <- function(df, col_select, name, sep=',', quote=FALSE, rowname=FALSE){
  write.table(df[col_select], name, sep=sep, quote=quote,row.names = rowname)
}

### read in file
dataDir ='/GCS_calling/GCS_merge_pseudo_strength.txt' # dataframe contains pseudo strengths, takes the union of all GCS
GCS_union_rep<-read.csv(file = file.path(getwd(), dataDir),header = TRUE) %>%
  as_tibble()

#prepare annotation for the heat map
sample_name = c('MOXI','LEVO','GEMI','CIP','NOR')
Drug <-factor(rep(sample_name,each=3),levels = sample_name)
meta <- data.frame(condition = colnames(GCS_union_rep)[2:16],Drug=Drug)

GCS_union_wide <- GCS_union_rep %>%
  melt(id.vars = 'pos',variable.name='condition',value.name='Strength')%>%
  left_join(meta)%>%
  group_by(pos,Drug) %>%
  summarise(avg_strength = mean(Strength))%>%
  spread('Drug','avg_strength')


# get the largest value column and row, sorted
MOXI_sort <-  GCS_union_wide %>%  
    filter( MOXI > LEVO & MOXI >GEMI & MOXI >CIP & MOXI >NOR)%>%  
    arrange(desc(MOXI))

LEVO_sort <-  GCS_union_wide %>%  
  filter( LEVO > MOXI & LEVO >GEMI & LEVO >CIP & LEVO >NOR)%>%  
  arrange(desc(LEVO))

GEMI_sort <-  GCS_union_wide %>%  
  filter( GEMI > MOXI & GEMI >LEVO &GEMI > CIP & GEMI >NOR)%>%  
  arrange(desc(GEMI))

CIP_sort <-  GCS_union_wide %>%  
  filter( CIP > MOXI & CIP > LEVO & CIP > GEMI & CIP > NOR)%>%  
  arrange(desc(CIP))

NOR_sort <-  GCS_union_wide %>%  
  filter( NOR > MOXI & NOR > LEVO & NOR > GEMI & NOR > CIP)%>%  
  arrange(desc(NOR))


save_df(LEVO_sort, c('pos','LEVO'), 'TopGCS_LEVO_sort.txt' ) #for motif plot
save_df(MOXI_sort, c('pos','MOXI'), 'TopGCS_MOXI_sort.txt' )
save_df(GEMI_sort, c('pos','GEMI'), 'TopGCS_GEMI_sort.txt' )
save_df(CIP_sort, c('pos','CIP'), 'TopGCS_CIP_sort.txt' )
save_df(NOR_sort, c('pos','NOR'), 'TopGCS_NOR_sort.txt' )



# heatmap -----------------------------------------------------------------
sort_bind <- rbind(MOXI_sort,LEVO_sort ,GEMI_sort,CIP_sort,NOR_sort) %>%
  column_to_rownames('pos')

### annotation color
annoCol <- brewer.pal(9, "Set1")[1:5]
names(annoCol) <- unique(meta$Drug)
annoCol <- list(Drug=annoCol)
heat_colors <- brewer.pal(9, "YlOrRd")

heatmap(as.matrix(sort_bind),Colv = NA, Rowv = NA)




# bar plot top GCS percentage -------------------------------------------------
data <- data.frame(Drug = factor(sample_name, levels = c('MOXI','LEVO','GEMI','CIP','NOR')), topGCS_ptg = c(11975/23272, 5493/23272, 4242/23272, 1375/23272,187/23272),
                   topGCS_ct = c(11975, 5493, 4242, 1375,187))

ap1 = ggplot(data , aes(x = Drug, y = topGCS_ptg))+
  geom_bar(stat = 'identity') +
  theme_classic()+
  ylab('% of Highest Cleavage Strength')+
  xlab('')


# number of top GCSs as genuine GCS to the specific FQ --------------------
###### input --------------------------------------------------------------
dataDir <- '/GCS_calling/GCS_extended_version/'
dfFiles <-  list.files(file.path(getwd(), dataDir), pattern="*.txt", full.name=T)
names(dfFiles) <- gsub(".txt","",basename(dfFiles))
df_list <- lapply(dfFiles, read.table, sep = ',',header = TRUE) #list of GCSs from 5 FQs

sort_list <- list(CIP_sort$pos,GEMI_sort$pos,LEVO_sort$pos,MOXI_sort$pos,NOR_sort$pos)
names(sort_list) <- c('CIP','GEMI','LEVO','MOXI','NOR')

hit_GCS <- c()
for (name in data$Drug){
  print(name)
  ct <- sum(sort_list[[name]] %in% df_list[[name]]$pos)
  print(ct)
  hit_GCS <- append(hit_GCS, ct)
}

data %>%
  mutate(real_GCS=hit_GCS, real_ptg = hit_GCS/topGCS_ct)%>%
  ggplot()+
  geom_bar(aes(x = Drug, y =real_ptg),stat = 'identity') +
  theme_classic()+
  ylab(bquote(atop('% Top Cleavage Strenth', 'Sites as GCSs')))+
  xlab('')


# top hits for motif construction -----------------------------------------

Top_Hits <- function(df_top,name){
  hit = df_top[c('pos',name)] %>% 
    filter(pos %in% df_list[[name]]$pos)
  return(hit)
}


save_df(Top_Hits(LEVO_sort,'LEVO'), c('pos','LEVO'), 'LEVO_sort_hit.txt') #for motif plot
save_df(Top_Hits(MOXI_sort,'MOXI'), c('pos','MOXI'), 'MOXI_sort_hit.txt')
save_df(Top_Hits(GEMI_sort,'GEMI'), c('pos','GEMI'), 'GEMI_sort_hit.txt')
save_df(Top_Hits(CIP_sort,'CIP'), c('pos','CIP'), 'CIP_sort_hit.txt')
save_df(Top_Hits(NOR_sort,'NOR'), c('pos','NOR'), 'NOR_sort_hit.txt')

