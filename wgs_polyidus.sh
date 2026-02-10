#!/bin/bash
python ~/tools/polyidus/src/polyidus.py -h

pip install -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple pandas
pip install -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple psutil
pip install -i https://mirrors.tuna.tsinghua.edu.cn/pypi/web/simple pysam

#/home/chengennong/tools/bowtie2/bowtie2-build -f genome.fa hg38
/home/chengennong/tools/bowtie2/bowtie2-build -f sb_aavr_car.fa sb

for sample in LV-1 LV-2 LV-3
do
cd /mnt/cgn/data/aavr/wgs_240902/$sample
mkdir -p polyidus && cd polyidus
python ~/tools/polyidus/src/polyidus.py \
    /mnt/cgn/data/ref/wgs_ref/hg38 \
    /mnt/cgn/data/ref/wgs_ref/lenti \
    --aligner bowtie2 \
    --fastq /mnt/cgn/data/aavr/wgs_240902/$sample/R1_qc.fastq.gz \
    /mnt/cgn/data/aavr/wgs_240902/$sample/R2_qc.fastq.gz \
    --outdir /mnt/cgn/data/aavr/wgs_240902/$sample/polyidus \
    --skip-alignment --virname LV
awk 'NR > 1' results/exactLVIntegrations.bed > results/exactLVIntegrations1.bed
done

for sample in AV-1 AV-2 AV-3
do
cd /mnt/cgn/data/aavr/wgs_240902/$sample
mkdir -p polyidus && cd polyidus
python ~/tools/polyidus/src/polyidus.py \
    /mnt/cgn/data/ref/wgs_ref/hg38 \
    /mnt/cgn/data/ref/wgs_ref/sb \
    --aligner bowtie2 \
    --fastq /mnt/cgn/data/aavr/wgs_240902/$sample/R1_qc.fastq.gz \
    /mnt/cgn/data/aavr/wgs_240902/$sample/R2_qc.fastq.gz \
    --outdir /mnt/cgn/data/aavr/wgs_240902/$sample/polyidus \
    --skip-alignment --virname SB
awk 'NR > 1' results/exactSBIntegrations.bed > results/exactSBIntegrations1.bed
done

#sequence
#CAGTTGAAGTCGGAAGTTTACATACACTTAAGTTG #sb_35
#CAGTTGAAGTCGGAAGTTTACATACACTTA #sb_30
#CAGTTGAAGTCGGAAGTTTACATACAC #e_short
#CAGTTGAAGTCGGAAGTTTACATACACCTTAGCCAAA #e

#ctgtactgggtctctctggttagaccagatctgag #lv_35

cat /mnt/cgn/data/aavr/wgs_240902/AV-1/polyidus/results/exactSBIntegrations.tsv
cat /mnt/cgn/data/aavr/wgs_240902/AV-1/polyidus/viral/ViralAligned_1.fastq |grep TATACAGTTGAAGTCG -B1
cat /mnt/cgn/data/aavr/wgs_240902/AV-1/polyidus/viral/ViralAligned_2.fastq |grep TATACAGTTGAAGTCG -B1

cat /mnt/cgn/data/ref/wgs_ref/lenti_aavr_car.fa |grep -i CTGTACTGGGTCTCTCTGGTTAGACCAGATCTGAG 
cat /mnt/cgn/data/aavr/wgs_240902/LV-*/polyidus/viral/ViralAligned_1.fastq | \
    grep -i ccttttagtcagtgtggaaaatctctagca | grep -v -i GTGGAAAATCTCTAGCAGTGGCGCCCG #| wc -l 3'LTR

cat /mnt/cgn/data/aavr/wgs_240902/LV-*/polyidus/viral/ViralAligned_2.fastq | \
    grep -i ctgtactgggtctctctggt | grep -v -i TCCGGACTGTACTGGGTCTCTCTGGT #| wc -l  5'LTR


cat /mnt/cgn/data/aavr/wgs_240902/LV-1/polyidus/results/LVIntegrationInfo.tsv
cat /mnt/cgn/data/aavr/wgs_240902/LV-1/polyidus/viral/ViralAligned_2.fastq |grep LH00169:446:22FKKWLT4:8:1208:14886:22264 -a1
#ctggaagggctaattcactcccaacgaagacaagatatccttgatctg
samtools faidx genome.fa chr1:7920820-7920836
