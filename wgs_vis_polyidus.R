#!/home/chengennong/tools/mambaforge/envs/genome/bin/R

library(tidyverse)
library(GenomicRanges)
library(genomation)
library(ggpubr)

#insertion read
wgs_res_csv_files <- list.files("/mnt/cgn/data/aavr/wgs_240902/",recursive = T,pattern = "Info.tsv",full.names = T)
wgs_res_csv_files <- wgs_res_csv_files[1:3]
random_res_csv_files <- read.table("/home/chengennong/ylp/others/lyj/aavr/result/Random.txt", header=T)
wgs_res_df <- (map_dfr(wgs_res_csv_files,read.table, sep="\t", header=T)) %>% 
  mutate(sample=str_extract(ViralFile,".V-."))

(wgs_res_df %>% filter(ChromViral=="Lenti") %>% select(PositionViral,ChromViral))
(wgs_res_df %>% filter(ChromViral=="SB") %>% select(PositionViral,ChromViral))
wgs_res_df$ChromViral
anno_dir <- "/home/chengennong/ylp/others/lyj/aavr/gene_coor"
features_bed_file = paste0(anno_dir,"/", (list.files(anno_dir)))
df_insertion_annot <- data.frame()

for (sample_i in unique(wgs_res_df$sample)){
  wgs_pd_res <- wgs_res_df %>% filter(sample==sample_i) %>% select(ChromHost,PositionHost,StrandHost, sample, Score)  %>%
    separate_rows(ends_with("Host"),sep = ", ") %>% 
    mutate(start=as.integer(PositionHost)-1,end=as.integer(PositionHost)) %>%
    mutate(strand=if_else(StrandHost=="Positive","+","-")) %>%
    mutate(readnames=rownames(.)) %>% mutate(seqname=ChromHost) %>%
    select(seqname,start, end, readnames, sample, strand)
    insertion_bed <- makeGRangesFromDataFrame(wgs_pd_res)
for (i in 1:length(features_bed_file)){
    features_bed <- readBed(features_bed_file[i],track.line=FALSE,remove.unusual=FALSE, zero.based=TRUE)
    insertion_annot_i <- annotateWithFeature(target=insertion_bed, feature=features_bed, intersect.chr = TRUE)
    features_name_v <- (features_bed_file[i] %>% str_split("/"))[[1]]; features_name <- features_name_v[length(features_name_v)] %>% str_remove(".bed")
    sample_name <- sample_i
    df_insertion_annot_i <- data.frame(number=insertion_annot_i@num.annotation, percent=insertion_annot_i@annotation) %>% 
      mutate(feature=features_name, sample=sample_name, total_number=sum(insertion_annot_i@num.annotation)) %>%
      filter(rownames(.) != "other")
    df_insertion_annot  <- rbind(df_insertion_annot,df_insertion_annot_i)
}}
wgs_res_df <- rbind(df_insertion_annot,random_res_csv_files)  %>% 
  separate(col = sample, into = c("type","id"), sep = "-",remove =  F) %>% 
  mutate(feature=str_replace(feature,"not","not "))

p_genic_wgs <- ggplot(wgs_res_df %>% filter(feature != "SafeHarbors_hg38_CLASSIC" & total_number != 1e5),
                      aes(x=feature, y=sample, fill=percent)) +
  geom_tile(color="grey")+ 
  geom_text(aes(label=round(percent,0)))+
  scale_fill_gradient(low = "brown", high = "white")+
  #scale_color_manual(values =c(rep("white",9)))+
  labs(title="Insertion read in genic features: WGS_PD")+theme_bw()+
  coord_flip()+
  rotate_x_text(angle = 45)
#ggsave(p_genic_wgs, filename="/home/chengennong/ylp/others/lyj/aavr/fig/p_genic_wgs_pd.pdf", width=5,height=3)

p_safe_wgs <- ggplot(wgs_res_df %>% filter(feature == "SafeHarbors_hg38_CLASSIC" & total_number != 1e5),
                  aes(x=sample, y=percent, fill=type)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  scale_fill_manual(values =c("#006400", "orange",  "grey"))+
  labs(title="Insertion read in safe harbors: WGS_PD")+theme_bw()+
  rotate_x_text(angle = 45)
#ggsave(p_safe_wgs, filename="/home/chengennong/ylp/others/lyj/aavr/fig/p_safe_wgs_pd.pdf", width=4.5,height=3)

stat_df <- wgs_res_df %>% filter(type !="Random") %>% select(sample,total_number, type) %>% unique
p_stat_PD <- ggplot(stat_df,
                  aes(x=sample, y=total_number, fill=type)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  scale_fill_manual(values =c("#006400", "orange",  "grey"))+
  labs(title="Insertion read number: WGS_PD")+theme_bw()+
  rotate_x_text(angle = 45)
#ggsave(p_stat_PD, filename="/home/chengennong/ylp/others/lyj/aavr/fig/p_stat_wgs_pd.pdf", width=4.5,height=3)
p_read <- ggarrange(p_stat_PD,p_genic_wgs,p_safe_wgs, ncol =1)


#insertion site
wgs_res_csv_files <- list.files("/mnt/cgn/data/aavr/wgs_240902/",recursive = T,pattern = "Info.tsv",full.names = T)
wgs_res_csv_files <- wgs_res_csv_files[1:3]
random_res_csv_files <- read.table("/home/chengennong/ylp/others/lyj/aavr/result/Random.txt", header=T)
wgs_res_df <- (map_dfr(wgs_res_csv_files,read.table, sep="\t", header=T)) %>% 
  mutate(sample=str_extract(ViralFile,".V-."))

(wgs_res_df %>% filter(ChromViral=="Lenti") %>% select(PositionViral,ChromViral))
(wgs_res_df %>% filter(ChromViral=="SB") %>% select(PositionViral,ChromViral))
wgs_res_df$ChromViral
anno_dir <- "/home/chengennong/ylp/others/lyj/aavr/gene_coor"
features_bed_file = paste0(anno_dir,"/", (list.files(anno_dir)))
df_insertion_annot <- data.frame()

for (sample_i in unique(wgs_res_df$sample)){
  wgs_pd_res <- wgs_res_df %>% filter(sample==sample_i) %>% select(ChromHost,PositionHost,StrandHost, sample, Score)  %>%
    separate_rows(ends_with("Host"),sep = ", ") %>% unique %>%
    mutate(start=as.integer(PositionHost)-1,end=as.integer(PositionHost)) %>%
    mutate(strand=if_else(StrandHost=="Positive","+","-")) %>%
    mutate(readnames=rownames(.)) %>% mutate(seqname=ChromHost) %>%
    select(seqname,start, end, readnames, sample, strand)
    insertion_bed <- makeGRangesFromDataFrame(wgs_pd_res)
for (i in 1:length(features_bed_file)){
    features_bed <- readBed(features_bed_file[i],track.line=FALSE,remove.unusual=FALSE, zero.based=TRUE)
    insertion_annot_i <- annotateWithFeature(target=insertion_bed, feature=features_bed, intersect.chr = TRUE)
    features_name_v <- (features_bed_file[i] %>% str_split("/"))[[1]]; features_name <- features_name_v[length(features_name_v)] %>% str_remove(".bed")
    sample_name <- sample_i
    df_insertion_annot_i <- data.frame(number=insertion_annot_i@num.annotation, percent=insertion_annot_i@annotation) %>% 
      mutate(feature=features_name, sample=sample_name, total_number=sum(insertion_annot_i@num.annotation)) %>%
      filter(rownames(.) != "other")
    df_insertion_annot  <- rbind(df_insertion_annot,df_insertion_annot_i)
}}
wgs_res_df <- rbind(df_insertion_annot,random_res_csv_files)  %>% 
  separate(col = sample, into = c("type","id"), sep = "-",remove =  F) %>% 
  mutate(feature=str_replace(feature,"not","not "))

p_genic_wgs_site <- ggplot(wgs_res_df %>% filter(feature != "SafeHarbors_hg38_CLASSIC" & total_number != 1e5),
                      aes(x=feature, y=sample, fill=percent)) +
  geom_tile(color="grey")+ 
  geom_text(aes(label=round(percent,0)))+
  scale_fill_gradient(low = "brown", high = "white")+
  #scale_color_manual(values =c(rep("white",9)))+
  labs(title="Insertion site in genic features: WGS_PD")+theme_bw()+
  coord_flip()+
  rotate_x_text(angle = 45)
#ggsave(p_genic_wgs, filename="/home/chengennong/ylp/others/lyj/aavr/fig/p_genic_wgs_pd.pdf", width=5,height=3)

p_safe_wgs_site <- ggplot(wgs_res_df %>% filter(feature == "SafeHarbors_hg38_CLASSIC" & total_number != 1e5),
                  aes(x=sample, y=percent, fill=type)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  scale_fill_manual(values =c("#006400", "orange",  "grey"))+
  labs(title="Insertion site in safe harbors: WGS_PD")+theme_bw()+
  rotate_x_text(angle = 45)
#ggsave(p_safe_wgs, filename="/home/chengennong/ylp/others/lyj/aavr/fig/p_safe_wgs_pd.pdf", width=4.5,height=3)

stat_df <- wgs_res_df %>% filter(type !="Random") %>% select(sample,total_number, type) %>% unique
p_stat_PD_site <- ggplot(stat_df,
                  aes(x=sample, y=total_number, fill=type)) +
  geom_bar(stat="identity", position=position_dodge())+ 
  scale_fill_manual(values =c("#006400", "orange",  "grey"))+
  labs(title="Insertion site number: WGS_PD")+theme_bw()+
  rotate_x_text(angle = 45)
#ggsave(p_stat_PD, filename="/home/chengennong/ylp/others/lyj/aavr/fig/p_stat_wgs_pd.pdf", width=4.5,height=3)
p_site <- ggarrange(p_stat_PD_site,p_genic_wgs_site,p_safe_wgs_site, ncol =1)

ggsave(ggarrange(p_read,p_site,ncol=2), filename="/home/chengennong/ylp/others/lyj/aavr/fig/p_pd_read_sbonly.pdf", width=8,height=8)
