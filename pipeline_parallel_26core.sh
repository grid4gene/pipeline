#!/bin/bash

set -x
ulimit -n 65535

#input file
export base='/home'

#input file
sp_dir='/mnt/disk_nvme/20190826A1B1-1'
#sample="S065_jiade-A_HGC-1911366_AHMYKWDSXX_S1_L002"
sample="S066_jiade-A_HGC-1911367_AHMYKWDSXX_S2_L002"

mkdir -p "${sp_dir}/input"
mkdir -p "${sp_dir}/output"

export R1="${sp_dir}/input/${sample}_R1_001.fastq.gz"
export R2="${sp_dir}/input/${sample}_R2_001.fastq.gz"

#reference file
export reference="$base/reference"
export hg19="$base/reference/hg19.fasta"
export dir_vcf="$base/reference"
export threads=26
export cores=26
export sp_name="gatk-3.8-fast-$sample"
export dir_analyzed_sample="${sp_dir}/output"
export dir_bin='/home/test/WGS_pipeline/TOOLS/bin'
export bed="$base/reference/Agilent_S06588914_Covered.bed"

# Align the reads
time bash -c " ${dir_bin}/bwa mem  -t $threads  $hg19  $R1 $R2 | ${dir_bin}/samtools view -@ $threads -S -b > $dir_analyzed_sample/$sp_name-bwa.bam"
# filter reads un-mapped reads and sort bam file
time bash -c "${dir_bin}/samtools view -@ $cores -b -h  -F 4 $dir_analyzed_sample/$sp_name-bwa.bam  > $dir_analyzed_sample/$sp_name-ali.bam"
# Sorting BAM
time bash -c " ${dir_bin}/samtools sort -m 8GB -@ $cores $dir_analyzed_sample/$sp_name-ali.bam  -o $dir_analyzed_sample/$sp_name-ali-sorted.bam"

# Add ReadGroup
time bash -c " java -Dsamjdk.compression_level=1 -jar ${dir_bin}/picard.jar AddOrReplaceReadGroups  I=$dir_analyzed_sample/$sp_name-ali-sorted.bam  O=$dir_analyzed_sample/$sp_name-ali-sorted-RG.bam  RGSM=$sp_name RGLB=project RGPL=illumina RGID=none RGPU=none VALIDATION_STRINGENCY=LENIENT"

# Remove duplicates
time bash -c " java -Dsamjdk.compression_level=1 -XX:+UseParallelGC -XX:ParallelGCThreads=2 -jar ${dir_bin}/picard.jar MarkDuplicates  REMOVE_DUPLICATES=FALSE CREATE_INDEX=TRUE VALIDATION_STRINGENCY=LENIENT MAX_RECORDS_IN_RAM=4000000 ASSUME_SORTED=TRUE I=$dir_analyzed_sample/$sp_name-ali-sorted-RG.bam  O=$dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam  METRICS_FILE=$dir_analyzed_sample/$sp_name-pacard.metrics"

# index BAM
time bash -c " ${dir_bin}/samtools index -@ $cores $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam"

# RealignerTargetCreator
time bash -c " java -Dsamjdk.compression_level=1 -jar $dir_bin/GenomeAnalysisTK.jar -nt $cores -T RealignerTargetCreator -known $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf  -known $dir_vcf/1000G_phase1.indels.hg19.sites.vcf  -R $hg19 -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam -o $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam.list"

# IndelRealigner
time bash -c " java -Dsamjdk.compression_level=1 -jar $dir_bin/GenomeAnalysisTK.jar  -T IndelRealigner -known $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known $dir_vcf/1000G_phase1.indels.hg19.sites.vcf -R $hg19 -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam -targetIntervals $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam.list -o $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.bam"


#BaseRecalibrator
#time bash -c " java -Dsamjdk.compression_level=1 -jar $dir_bin/GenomeAnalysisTK.jar  -nct $cores -T BaseRecalibrator -l INFO  -R $hg19   -knownSites $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $dir_vcf/1000G_phase1.indels.hg19.sites.vcf  -I  $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.bam -o  $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.grp"
# GATK PrintReads
#time bash -c " java -Dsamjdk.compression_level=1 -jar $dir_bin/GenomeAnalysisTK.jar -nct 8 -T PrintReads -R $hg19 -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.bam -BQSR $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.grp -o $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam"

#create script for Scatter Gather for bqsr and HaplotypeCaller
time python create_script.py
#BaseRecalibrator
time bash -c "parallel -j6 < ./script/bqsr.sh"
time ./script/merge_bqsr.sh
# GATK PrintReads
time bash -c "parallel -j6 < ./script/apply.sh"
time ./script/merge_bam.sh

# Filter mapping quality <10
#time bash -c "${dir_bin}/samtools view -@ $threads -h -q 10 $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal-0.bam |${dir_bin}/samtools view  -@ $threads -h -Sb > $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam"

# index BAM
time bash -c " ${dir_bin}/samtools index -@ $cores $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam"

# GATK UnifiedGenotyper
time bash -c " java -Dsamjdk.compression_level=1 -jar $dir_bin/GenomeAnalysisTK.jar -nt $cores -T UnifiedGenotyper -R $hg19  -L $bed  -metrics $dir_analyzed_sample/$sp_name-snps.metrics -stand_call_conf 10.0 -dcov 20000 -glm BOTH -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o $dir_analyzed_sample/$sp_name-Unified-SNP-INDLE.vcf"
#time java -Dsamjdk.compression_level=1 -jar $dir_bin/GenomeAnalysisTK.jar -nt $threads -T UnifiedGenotyper -R $hg19  -metrics $dir_analyzed_sample/$sp_name-snps.metrics -stand_call_conf 10.0 -stand_emit_conf 5.0 -dcov 20000 -glm BOTH -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o $dir_analyzed_sample/$sp_name-Unified-SNP-INDLE.vcf

#GARK HaplotypeCaller
#time java -Dsamjdk.compression_level=1 -jar $dir_bin/GenomeAnalysisTK.jar -nct  $threads -T HaplotypeCaller -R $hg19 -L $bed  -stand_call_conf 10.0 -stand_emit_conf 5.0  -minPruning 3  -mbq 5  -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o $dir_analyzed_sample/$sp_name-Haploy-SNP-INDLE.vcf
#time bash -c "java -Dsamjdk.compression_level=1 -jar $dir_bin/GenomeAnalysisTK.jar -nct 4 -T HaplotypeCaller -R $hg19 -stand_call_conf 10.0 -minPruning 3  -mbq 5  -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o $dir_analyzed_sample/$sp_name-Haploy-SNP-INDLE.vcf"
time bash -c "parallel -j6 < ./script/htc.sh"
time ./script/merge_htc.sh
