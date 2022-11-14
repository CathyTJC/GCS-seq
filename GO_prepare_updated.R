###### Add GCS annotation/ prepare for functional enrichment analysis
source('utils_data.R')

write_geneID <- function(FQ_ano, thereshold,output_name){
  # write file for DAVID
  FQ_id<-FQ_ano %>% group_by(gene_name)%>%
    summarise(gene_strength = sum(avg_strength),gene_id=gene_id2)%>%
    unique()%>%
    filter(gene_strength>thereshold)
  write(FQ_id$gene_id, paste0(output_name, '.txt'))
  
  return (FQ_id)
}

# Input:Read in GCS (include Mu region) -------------------------------------------------------------------

# read in GCS file with genome cleavage strength
dataDir <- '/data/GCS/GCS_extended_version/'
dfFiles <-  list.files(file.path(getwd(), dataDir), pattern="*.txt", full.name=T)
names(dfFiles) <- gsub(".txt","",basename(dfFiles))
df_list <- lapply(dfFiles, read.table, sep = ',',header = TRUE)


#read in  annotation file: Mu-ori genome
.anoPath <- '/data/REFERENCE/GCF_000005845.2_ASM584v2_genomic_insert_hierexon.gffread.gtf'
.colselect = c("start", "end","gene_id","gene_name")
muori_ref<-readAnnotaion(AnnotationPath =.anoPath,col_to_select=.colselect)

# add annotation (within ORF) ---------------------------------------------
FQ_ano_ls <- lapply(df_list, add_ano, muori_ref)

#Write files---------------------------------------------
write_geneID(FQ_ano=FQ_ano_ls$CIP, thereshold=1,output_name='geneID_CIP1')
write_geneID(FQ_ano=FQ_ano_ls$MOXI, thereshold=1,output_name='geneID_MOXI1')
write_geneID(FQ_ano=FQ_ano_ls$LEVO, thereshold=1,output_name='geneID_LEVO1')
write_geneID(FQ_ano=FQ_ano_ls$NOR, thereshold=1,output_name='geneID_NOR1')
write_geneID(FQ_ano=FQ_ano_ls$GEMI, thereshold=1,output_name='geneID_GEMI1')

#all genes
write_geneID(FQ_ano=FQ_ano_ls$CIP, thereshold=0,output_name='geneID_CIP')
write_geneID(FQ_ano=FQ_ano_ls$MOXI, thereshold=0,output_name='geneID_MOXI')
write_geneID(FQ_ano=FQ_ano_ls$LEVO, thereshold=0,output_name='geneID_LEVO')
write_geneID(FQ_ano=FQ_ano_ls$NOR, thereshold=0,output_name='geneID_NOR')
write_geneID(FQ_ano=FQ_ano_ls$GEMI, thereshold=0,output_name='geneID_GEMI')

