################
# GCS hierarchical clustering heat map
###### 

library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ComplexHeatmap)

### read in file
dataDir ='data/GCS_calling/GCS_merge_pseudo_strength.txt' # dataframe contains pseudo strengths, takes the union of all GCS
GCS_union_rep<-read.csv(file = file.path(getwd(), dataDir),header = TRUE,row.names = 1)

dataDir ='data/GCS_calling/GCS_merge_shared.csv' 
GCS_shared <- read.csv(file = file.path(getwd(), dataDir),header = TRUE,row.names = 2) %>% select (-X)
colnames(GCS_shared) =colnames(GCS_union_rep)
  
#prepare annotation for the heat map
sample_name = c('MOXI','LEVO','GEMI','CIP','NOR')
Drug <-factor(rep(sample_name,each=3),levels = sample_name)
meta <- data.frame(row.names = colnames(GCS_union_rep),Drug=Drug)
### annotation color
annoCol <- brewer.pal(9, "Set1")[1:5]
names(annoCol) <- unique(meta$Drug)
annoCol <- list(Drug=annoCol)

### set manual color scale
heat_colors <- brewer.pal(9, "YlOrRd")
quantile_breaks <- function(xs, n = 20) {
  breaks <- quantile(xs, probs = seq(0, 1, length.out = n))
  breaks[!duplicated(breaks)]
}



### plot: union of GCSs
mat_breaks <- quantile_breaks(as.matrix(GCS_union_rep), n = 10)
ComplexHeatmap::pheatmap(GCS_union_rep, name = ' ',
         color = colorRampPalette(heat_colors)(length(mat_breaks)-1),cluster_rows =T,show_rownames = F,
         breaks = mat_breaks,border_color = NA,annotation_col =meta,heatmap_legend_param = list(at = round(mat_breaks,2), color_bar = "discrete"),
         annotation_colors = annoCol,main = "Quantile Color Scale") 
# pheatmap(GCS_union_rep,color=heat_colors,cluster_rows =T,show_rownames = F,scale = 'row', border_color = NA,annotation=meta)


# shared GCS heatmap ------------------------------------------------------
mat_breaks <- quantile_breaks(as.matrix(GCS_shared), n = 10)
ComplexHeatmap::pheatmap(GCS_shared, name=' ',
         color = colorRampPalette(heat_colors)(length(mat_breaks)-1),cluster_rows =T,show_rownames = F,
         breaks = mat_breaks,border_color = NA,annotation_col =meta,heatmap_legend_param = list(at = round(mat_breaks,2), color_bar = "discrete"),
         annotation_colors = annoCol, main = "956 Shared GCSs") 


