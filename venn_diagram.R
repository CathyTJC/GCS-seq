# Author: "Juechun Tang"
# Venn diagram of GCS from LEVO, MOXI, NOR, CIP, GEMI 30mins treatment

library(tidyverse)
library(VennDiagram)
library(RColorBrewer)
library(venn)
library(gridExtra)
library(ggpmisc)
library(ggpubr)

# read in GCSs
setwd("/Volumes/TJC/GCS_calling_files")
dataDir <- "HF_GCS_strength_0408/"
dfFiles <- dir(file.path(getwd(), dataDir), pattern="*.txt", full.name=T) #raw GCS files
names(dfFiles) <- gsub(".txt","",basename(dfFiles))
CIP <-read.table(file = dfFiles[1], sep = ',',header = TRUE)$pos
GEMI <-read.table(file = dfFiles[2], sep = ',',header = TRUE)$pos 
LEVO <-read.table(file = dfFiles[3], sep = ',',header = TRUE)$pos
MOXI <-read.table(file = dfFiles[4], sep = ',',header = TRUE)$pos 
NOR <-read.table(file = dfFiles[5], sep = ',',header = TRUE)$pos

# 5 fq comparison---------------------------------------------------
venn.result <- venn(list(LEVO =LEVO, MOXI=MOXI, NOR=NOR, CIP = CIP, GEMI= GEMI), ilabels = TRUE, 
                    zcolor = "style", opacity = 0.1)
