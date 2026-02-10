#!/home/chengennong/tools/mambaforge/envs/genome/bin/R

library(tidyverse)
library(circlize)
library(ggh4x)
library(Biostrings)
# 读取FASTA文件
#fasta_file <- "/mnt/cgn/data/ref/wgs_ref/genome.fa"
#seqs <- readDNAStringSet(fasta_file)
#filtered_seqs <- seqs[grep("^chr(?!M)", names(seqs), perl = TRUE)]
#ref <- rbind(data.frame(Genome=names(filtered_seqs),Length=width(filtered_seqs),spe="Human"))
    #data.frame(Genome="LV",Length=width(readDNAStringSet("/mnt/cgn/data/ref/wgs_ref/lenti_aavr_car.fa")),spe="LV"),
    #data.frame(Genome="SB",Length=width(readDNAStringSet("/mnt/cgn/data/ref/wgs_ref/sb_aavr_car.fa")),spe="SB"))

#sb
bed <- read_tsv("/mnt/cgn/data/aavr/wgs_240902/AV-1/polyidus/results/exactSBIntegrations.bed");bed
bed <- read_tsv("/mnt/cgn/data/aavr/wgs_240902/AV-1/polyidus/results/SBIntegrationInfo.tsv")
bed <- bed %>% select(ChromHost, PositionHost, StrandHost, ChromViral, NumberReads, PositionViral) %>% 
    mutate(Start_host=as.integer(str_split_fixed(PositionHost,pattern=", ", n=Inf)[,1]),End_host=as.integer(str_split_fixed(PositionHost,pattern=", ", n=Inf)[,1])+1) %>% 
    select(ChromHost,Start_host,End_host,ChromViral);bed

pdf("/home/chengennong/ylp/others/lyj/aavr/fig/cir_plot_sb.pdf", width=5, height=5)
circos.clear()  
circos.par(gap.after = c(rep(1,23),0), start.degree = 90)  
circos.initializeWithIdeogram(species = "hg38", plotType = NULL)  
#circos.initializeWithIdeogram(factors=ref$Genome,          
#                  xlim=matrix(c(rep(0,dim(ref)[1]),ref$Length),ncol=2))
#circos.genomicLabels(bed, labels.column = 4, side = "outside",
#                     connection_height = 0.1, labels.side = "clockwise")  
set_track_gap(mm_h(0.3))  
circos.track(ylim = c(0,1), track.height = 0.05, bg.border = "white", panel.fun = function(x, y) {
  sector.name <- CELL_META$sector.index  
  if(sector.name %in% paste0("chr", c(1:22, "X", "Y"))) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 1, sector.name, cex = 0.5,  
                facing = "inside", niceFacing = TRUE, adj = c(0.5, 0.5), col = "black") 
  }
})
circos.genomicIdeogram(track.height = mm_h(3))
df_repeated <- data.frame(ChromHost = "chr1", Start_host=0,
                          End_host=1, ChromViral="SB")[rep(1, dim(bed)[1]), ]
circos.genomicLink(df_repeated,col="#006400",
                   bed)
#circos.genomicLink(read_tsv("link1.tsv") %>% filter(col=="green"),
#                   read_tsv("link2.tsv") %>% filter(col=="green"),col="#0B775E")
dev.off()

#lv
bed <- read_tsv("/mnt/cgn/data/aavr/wgs_240902/LV-1/polyidus/results/exactLVIntegrations.bed");bed
bed <- read_tsv("/mnt/cgn/data/aavr/wgs_240902/LV-1/polyidus/results/LVIntegrationInfo.tsv")
bed <- bed %>% select(ChromHost, PositionHost, StrandHost, ChromViral, NumberReads, PositionViral) %>% 
    mutate(Start_host=as.integer(str_split_fixed(PositionHost,pattern=", ", n=Inf)[,1]),End_host=as.integer(str_split_fixed(PositionHost,pattern=", ", n=Inf)[,1])+1) %>% 
    select(ChromHost,Start_host,End_host,ChromViral);bed
pdf("/home/chengennong/ylp/others/lyj/aavr/fig/cir_plot_lv.pdf", width=5, height=5)
circos.clear() 
circos.par(gap.after = c(rep(1,23),0), start.degree = 90)  
circos.initializeWithIdeogram(species = "hg38", plotType = NULL)  

set_track_gap(mm_h(0.3)) 
circos.track(ylim = c(0,1), track.height = 0.05, bg.border = "white", panel.fun = function(x, y) {
  sector.name <- CELL_META$sector.index  
  if(sector.name %in% paste0("chr", c(1:22, "X", "Y"))) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 1, sector.name, cex = 0.5,  
                facing = "inside", niceFacing = TRUE, adj = c(0.5, 0.5), col = "black") 
  }
})
circos.genomicIdeogram(track.height = mm_h(3))
df_repeated <- data.frame(ChromHost = "chr1", Start_host=0,
                          End_host=1, ChromViral="Lenti")[rep(1, dim(bed)[1]), ]
circos.genomicLink(df_repeated,col="orange",
                   bed)
#circos.genomicLink(read_tsv("link1.tsv") %>% filter(col=="green"),
#                   read_tsv("link2.tsv") %>% filter(col=="green"),col="#0B775E")
dev.off()