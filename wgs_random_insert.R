
args <- commandArgs()
#print(args[6])
#print(args[7])
##print(args[c(8:length(args))])
#print(args[8])

input_bed_file = args[6]
#features_bed_file = args[c(8:length(args))]
anno_dir = args[7]
input_arm_file = args[8]
Arm_3 = args[9]
Arm_5 = args[10]
out_csv_file = args[11]

library(genomation)
library(GenomicRanges)
library(dplyr) 
library(stringr) 
packageVersion("GenomicRanges")


hg38_length_df <- read.table("/mnt/cgn/data/ref/wgs_ref/genome.fa.fai")
hg38_length <- regioneR::getGenomeAndMask(hg38_length_df %>% select(V1,V2) %>% filter(str_detect(V1,"chr")) %>% filter(!str_detect(V1,"chrM")), mask=NA)
anno_dir <- "/home/chengennong/ylp/others/lyj/aavr/gene_coor"
features_bed_file = paste0(anno_dir,"/", (list.files(anno_dir)))
df_insertion_annot <- data.frame()
for (nregion in c(1e4,1e5)){
for (seed_i in 1:3){
    set.seed(seed=seed_i)
    random_site <- regioneR::createRandomRegions(nregions=nregion, length.mean=2, length.sd=0, genome=hg38_length$genome, mask=NULL, non.overlapping=TRUE)
    #insertion_annot_r <- annotateWithFeature(target=random_site, feature=features_bed, intersect.chr = TRUE)
    #print(insertion_annot_r)
    bed_df <- as.data.frame(random_site);bed_df$name <- ".";bed_df$score <- "0"
    bed_df_bed <- bed_df %>% select(seqnames, start, end, name, score, strand)
    write.table(bed_df_bed, file = paste0("/home/chengennong/ylp/others/lyj/aavr/result/random",seed_i,".bed"), sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE)
for (i in 1:length(features_bed_file)){
    features_bed <- readBed(features_bed_file[i],track.line=FALSE,remove.unusual=FALSE, zero.based=TRUE)
    insertion_annot_i <- annotateWithFeature(target=random_site, feature=features_bed, intersect.chr = TRUE)
    features_name_v <- (features_bed_file[i] %>% str_split("/"))[[1]]; features_name <- features_name_v[length(features_name_v)] %>% str_remove(".bed")
    sample_name <- paste0("Random-",seed_i)
    df_insertion_annot_i <- data.frame(number=insertion_annot_i@num.annotation, percent=insertion_annot_i@annotation) %>% 
      mutate(feature=features_name, sample=sample_name, total_number=sum(insertion_annot_i@num.annotation)) %>%
      filter(rownames(.) != "other")
    df_insertion_annot  <- rbind(df_insertion_annot,df_insertion_annot_i)
}
}
}
df_insertion_annot

write.table(df_insertion_annot, "/home/chengennong/ylp/others/lyj/aavr/result/Random.txt", row.names= F)


