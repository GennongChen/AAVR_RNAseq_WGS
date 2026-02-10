#!/home/chengennong/tools/mambaforge/envs/atac/bin/R

library(GenomicFeatures)
#library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library()
library(ChIPseeker)
library(ggplot2)
gtf="/mnt/cgn/data/ref/human/gtf/genes.gtf"
genome="/mnt/cgn/data/ref/wgs_ref/genome.fa"
genecode.txdb<- txdbmaker::makeTxDbFromGFF(gtf, format="gtf", genome)

files_virus <- list.files("/mnt/cgn/data/aavr/wgs_240902",pattern="Integrations1.bed",recursive=T,full.names=T)
files_random <- list.files("/home/chengennong/ylp/others/lyj/aavr/result",pattern=".bed",recursive=T,full.names=T)
#files <- "/mnt/cgn/data/aavr/wgs_240902/AV-1/polyidus/results/exactSBIntegrations1.bed"
files <- c(files_virus, files_random)
print(files)
files <- files[c(1,2,3,7,8,9)]
options(ChIPseeker.ignore_1st_exon = TRUE)
options(ChIPseeker.ignore_1st_intron = TRUE)
options(ChIPseeker.ignore_downstream = TRUE)
#options(ChIPseeker.ignore_promoter_subcategory = TRUE)
peakAnnoList <- lapply(files, annotatePeak, tssRegion=c(-3000, 3000), 
                         TxDb=genecode.txdb, level = "transcript", annoDb="org.Hs.eg.db",
                         sameStrand = FALSE, ignoreOverlap = FALSE, 
                         ignoreDownstream = TRUE,
                         overlap = "all")

#peakAnno <- annotatePeak(files, tssRegion=c(-3000, 3000), 
#                         TxDb=genecode.txdb, level = "transcript", annoDb="org.Hs.eg.db",
#                         sameStrand = FALSE, ignoreOverlap = FALSE, #annotation = c("Promoter-TSS", "Exon", "Intron"),
#                         ignoreDownstream = TRUE, overlap = "all")
#names(peakAnnoList) <- c("AV-1","AV-2","AV-3","LV-1","LV-2","LV-3","Random-1","Random-2","Random-3")
#pdf("/home/chengennong/ylp/others/lyj/aavr/fig/p_anno.pdf", width=5,height=3)
names(peakAnnoList) <- c("AV-1","AV-2","AV-3","Random-1","Random-2","Random-3")
pdf("/home/chengennong/ylp/others/lyj/aavr/fig/p_anno_sbonly.pdf", width=5,height=2.5)
plotAnnoBar(peakAnnoList)
dev.off()

