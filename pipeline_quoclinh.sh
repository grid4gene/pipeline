#!/bin/bash

set -x
ulimit -n 65535

#default hdi
base='/home'
#number of threads
threads=26

#input file
sp_dir='/mnt/disk_nvme/20190826A1B1-1'
#sample="S065_jiade-A_HGC-1911366_AHMYKWDSXX_S1_L002"
sample="S066_jiade-A_HGC-1911367_AHMYKWDSXX_S2_L002"

R1="${sp_dir}/input/${sample}_R1_001.fastq.gz"
R2="${sp_dir}/input/${sample}_R2_001.fastq.gz"
#output
sp_name='gatk-3.8-origin-'$sample
dir_analyzed_sample="${sp_dir}/output"

#reference file
hg19="$base/reference/hg19.fasta"
dir_vcf="$base/reference"

#Tool dir
dir_bin='/home/test/WGS_pipeline/TOOLS/bin'
bed="$base/reference/Agilent_S06588914_Covered.bed"


time ${dir_bin}/bwa mem  -t $threads  $hg19  $R1 $R2 | ${dir_bin}/samtools view -S -b > ${dir_analyzed_sample}/${sp_name}-bwa.bam

# filter reads un-mapped reads and sort bam file
time ${dir_bin}/samtools view -h  -F 4 ${dir_analyzed_sample}/$sp_name-bwa.bam | ${dir_bin}/samtools view -Sb > ${dir_analyzed_sample}/$sp_name-ali.bam

# Sorting BAM
time ${dir_bin}/samtools sort ${dir_analyzed_sample}/$sp_name-ali.bam  -o ${dir_analyzed_sample}/$sp_name-ali-sorted.bam

# Add ReadGroup
time java -jar ${dir_bin}/picard.jar AddOrReplaceReadGroups  I=${dir_analyzed_sample}/$sp_name-ali-sorted.bam  O=${dir_analyzed_sample}/$sp_name-ali-sorted-RG.bam  RGSM=$sp_name RGLB=project RGPL=illumina RGID=none RGPU=none VALIDATION_STRINGENCY=LENIENT

# Remove duplicates
time java -jar ${dir_bin}/picard.jar MarkDuplicates  REMOVE_DUPLICATES=FALSE CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 ASSUME_SORTED=TRUE I=${dir_analyzed_sample}/$sp_name-ali-sorted-RG.bam  O=${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup.bam  METRICS_FILE=${dir_analyzed_sample}/$sp_name-pacard.metrics

# index BAM
time ${dir_bin}/samtools index ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup.bam

# RealignerTargetCreator
time java -jar $dir_bin/GenomeAnalysisTK.jar -nt $threads -T RealignerTargetCreator -known $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf  -known $dir_vcf/1000G_phase1.indels.hg19.sites.vcf  -R $hg19 -I ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup.bam -o ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup.bam.list

# IndelRealigner
time java -jar $dir_bin/GenomeAnalysisTK.jar  -T IndelRealigner -known $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known $dir_vcf/1000G_phase1.indels.hg19.sites.vcf -R $hg19 -I ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup.bam -targetIntervals ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup.bam.list -o ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned.bam

#BaseRecalibrator
time java -jar $dir_bin/GenomeAnalysisTK.jar  -nct $threads -T BaseRecalibrator -l INFO  -R $hg19   -knownSites $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $dir_vcf/1000G_phase1.indels.hg19.sites.vcf  -I  ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned.bam -o  ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned.grp

# GATK PrintReads
time java -jar $dir_bin/GenomeAnalysisTK.jar -nct $threads -T PrintReads -R $hg19 -I ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned.bam -BQSR ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned.grp -o ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam

# Filter mapping quality <10
#time ${dir_bin}/samtools view -h  -q 1 ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned-recal-0.bam |samtools view -h -Sb > ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam

# index BAM
time ${dir_bin}/samtools index ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam

# GATK UnifiedGenotyper
time java -jar $dir_bin/GenomeAnalysisTK.jar -nt $threads -T UnifiedGenotyper -R $hg19  -L $bed  -metrics ${dir_analyzed_sample}/$sp_name-snps.metrics -stand_call_conf 10.0 -dcov 20000 -glm BOTH -I ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o ${dir_analyzed_sample}/$sp_name-Unified-SNP-INDLE.vcf
#time java -jar $dir_bin/GenomeAnalysisTK.jar -nt $threads -T UnifiedGenotyper -R $hg19  -metrics ${dir_analyzed_sample}/$sp_name-snps.metrics -stand_call_conf 10.0 -stand_emit_conf 5.0 -dcov 20000 -glm BOTH -I ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o ${dir_analyzed_sample}/$sp_name-Unified-SNP-INDLE.vcf

#GARK HaplotypeCaller
#time java -jar $dir_bin/GenomeAnalysisTK.jar -nct  $threads -T HaplotypeCaller -R $hg19 -L $bed  -stand_call_conf 10.0 -stand_emit_conf 5.0  -minPruning 3  -mbq 5  -I ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o ${dir_analyzed_sample}/$sp_name-Haploy-SNP-INDLE.vcf
time java -jar $dir_bin/GenomeAnalysisTK.jar -nct  $threads -T HaplotypeCaller -R $hg19 -stand_call_conf 10.0 -minPruning 3  -mbq 5  -I ${dir_analyzed_sample}/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o ${dir_analyzed_sample}/$sp_name-Haploy-SNP-INDLE.vcf
