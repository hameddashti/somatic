source("https://bioconductor.org/biocLite.R")
biocLite("SomaticSignatures")

library(BSgenome.Hsapiens.UCSC.hg19)
library(SomaticSignatures)
library(GenomicRanges)
library(readr)
library(ggplot2)
library(LearnBayes)
library(dplyr)  
library(reshape)
library(VariantAnnotation)
library(data.table)
library(parallel)
library(tidyr)
library(csv)


Cancer<- data.frame(read_tsv("C://Users/my asus/Desktop/simple_somatic_mutation.open.tsv",col_names = TRUE))
gene_list <- suppressWarnings(read_excel('../../Data/genomic_annotation/FANTOM5_gene_list.xlsx'))
annot <- as.data.frame(gene_list[,c(1,3,4,5,6,8)])
colnames(annot) <- c('id','chr','start','end','strand','geneClass')
tmp<- as.data.frame(Cancer[,c(1,5,9,10,11,12,16,17)]) 
tmp1<- data.frame(subset(tmp, tmp$icgc_mutation_id !="NA" ))
tmp3<- data.frame(subset(tmp1, nchar(tmp1$mutated_from_allele)==1 & nchar(tmp1$mutated_to_allele)==1
                         & tmp1$mutated_from_allele !='-'& tmp1$mutated_to_allele != "-"
                         & tmp1$mutated_from_allele !='_'& tmp1$mutated_to_allele != "_"
                         & tmp1$chromosome != "MT" & tmp1$chromosome != "M") ) 

tmp2<- tmp3[!duplicated(tmp3),]

tmp2$chromosome_strand[which(tmp2$chromosome_strand == '1')] <- '+'
tmp2$chromosome_strand[which(tmp2$chromosome_strand == '2')] <- '-'

tmp2$chromosome <- paste0("chr",tmp2$chromosome)

gr <- makeGRangesFromDataFrame(tmp2,keep.extra.columns = T,
                               seqnames.field = 'chromosome',
                               start.field='chromosome_start',
                               end.field = 'chromosome_end',
                               strand.field = 'chromosome_strand')


idx <- match(c('mutated_from_allele','mutated_to_allele','icgc_sample_id'),names(mcols(gr)))
mcols(gr) <- cbind(mcols(gr)[idx],mcols(gr)[-idx])

  vr <- makeVRangesFromGRanges(gr,ref.field='mutated_from_allele',
                               alt.field='mutated_to_allele',
                               sampleNames.field = 'icgc_sample_id',
                               keep.extra.columns = T)
  
  vr <- mutationContext(vr,Hsapiens)
  
  
  variations <- data.frame(icgc_samlpe_id = mcols(gr)$icgc_sample_id,
                           chromosome = as.character(seqnames(vr)),
                           position = start(vr),
                           strand = as.character(strand(gr)),
                           motif = paste0(as.character(mcols(vr)$alteration),
                                               '-',as.character(mcols(vr)$context)))
  
  
  x<- as.data.frame(variations[,c(1,5,2,3)],header=T)
  
  df1 %>% inner_join(annot, "chr") %>% 
  mutate(geneID_motif = paste(id, motif, sep = ","),
         n = if_else(position >= start & position <= end, 1, 0)) %>% 
  select(icgc_samlpe_id, geneID_motif, n) %>%
  group_by(icgc_samlpe_id, geneID_motif) %>% 
  summarise(n = sum(n)) %>%
  spread(key = geneID_motif, value = n, fill = 0)

  
  
  

  write.csv(df1, "C://Users/my asus/Desktop/mydata.csv")  
  
  
 




