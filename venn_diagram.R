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


# 4Mu-LEVO vs LEVO --------------------------------------------------------
LEVO4Mu <-read.table("HF_GCS_strength_0408/4-Mu/4Mu-LEVO-strength.txt", sep = ',',header = TRUE)$mod_GCS

# Moxi and Levo Venn diagram
myCol <- brewer.pal(4, "Pastel2")
# Chart
venn.diagram(
  x = list(LEVO, LEVO4Mu),
  category.names = c("LEVO" , "LEVO-4Mu "),
  filename = '4Mu-LEVO_overlap.png',
  print.mode=c("raw","percent"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)


# compare each replicate --------------------------------------------------
setwd("/Volumes/TJC/GCS_calling_files/LEVO_compare_strength/")
LEVO <-read.table('LEVO.txt', sep = ',',header = TRUE)%>%
  select('pos','strength_rep1','strength_rep2','strength_rep3')
LEVO4Mu <- read.table("4Mu-LEVO-corrected.txt",sep = ',', header = TRUE)%>%
  select('mod_GCS','strength_rep1','strength_rep2','strength_rep3')
  
LEVOrep1 <-LEVO$pos[which(!is.na(LEVO$strength_rep1))]
LEVOrep2 <-LEVO$pos[which(!is.na(LEVO$strength_rep2))]
LEVOrep3 <-LEVO$pos[which(!is.na(LEVO$strength_rep3))]

LEVO4Murep1 <-LEVO4Mu$mod_GCS[which(!is.na(LEVO4Mu$strength_rep1))]
LEVO4Murep2 <-LEVO4Mu$mod_GCS[which(!is.na(LEVO4Mu$strength_rep2))]
LEVO4Murep3 <-LEVO4Mu$mod_GCS[which(!is.na(LEVO4Mu$strength_rep3))]

myCol <- brewer.pal(4, "Pastel2")
# Chart

s1 = venn.diagram(
  x = list(LEVOrep1, LEVOrep2),
  category.names = c("LEVOrep1" , "LEVO-rep2 "),
  filename = 'txt.png',
  print.mode=c("raw"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)


s2=venn.diagram(
  x = list(LEVOrep1, LEVOrep3),
  category.names = c("LEVOrep1" , "LEVO-rep3 "),
  filename = 'l13.png',
  print.mode=c("raw"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

s3 = venn.diagram(
  x = list(LEVOrep2, LEVOrep3),
  category.names = c("LEVOrep2" , "LEVO-rep3 "),
  filename = 'l23.png',
  print.mode=c("raw"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

s4 = venn.diagram(
  x = list(LEVO4Murep1, LEVO4Murep2),
  category.names = c("LEVO4Murep1" , "LEVO4Murep2 "),
  filename = '4l12.png',
  print.mode=c("raw"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

s5 = venn.diagram(
  x = list(LEVO4Murep1, LEVO4Murep3),
  category.names = c("LEVO4Murep1" , "LEVO4Murep3 "),
  filename = '4l13.png',
  print.mode=c("raw"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

s6 = venn.diagram(
  x = list(LEVO4Murep2, LEVO4Murep3),
  category.names = c("LEVO4Murep2" , "LEVO4Murep3 "),
  filename = '4l23.png',
  print.mode=c("raw"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

s7 = venn.diagram(
  x = list(LEVOrep1, LEVO4Murep1),
  category.names = c("LEVOrep1" , "LEVO4Murep1 "),
  filename = 's7.png',
  print.mode=c("raw"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

s15= venn.diagram(
  x = list(LEVOrep3, LEVO4Murep3),
  category.names = c("LEVOrep3" , "LEVO4Murep3 "),
  filename = 's15.png',
  print.mode=c("raw"),
  output=TRUE,
  
  # Output features
  imagetype="tiff" ,
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = myCol[2:3],
  
  # Numbers
  fontface = "bold",
  # Set names
  cat.fontface = "bold",
  cat.default.pos = "outer"
)

s_sum<-ggarrange(s1,s2,s3,s4,s5,s6,s7,s8,s9,s10,s11,s12,s13,s14,s15, ncol=3, nrow=5)
ggsave("venn.png", units="in", s_sum, width=7, height=8, dpi=1200)
