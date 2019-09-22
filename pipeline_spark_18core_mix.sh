#!/bin/bash
set -x

export base='/home'
export input='/home/input'
export output='/home/output'
export reference='/home/reference'
export tools='/home/tools'


#input file
#export sample="412"
#export R1="$input/${sample}_1.fastq.gz"
#export R2="$input/${sample}_2.fastq.gz"

export sample="NA12878"
export R1="$input/${sample}_R1_001.fastq.gz"
export R2="$input/${sample}_R2_001.fastq.gz"


#reference file
export hg19="$reference/hg19.fasta"
export dir_vcf="$reference"
export bed="$reference/Agilent_S06588914_Covered.bed"
#output file
export dir_analyzed_sample="$output"
export sp_name="pipe-3.8-fast-"$sample
#tools
export dir_bin="$tools"

export threads=18
export cores=18

ReadGroup="@RG\tID:${sample}_${cores}T\tLB:${sample}\tSM:${sample}\tPL:ILLUMINA"

export hdfs_url=localhost:9000/user/root
export spark_tools="/mnt/disk_a/spark_tools"

export spark_cmd="-- --spark-runner SPARK --spark-master yarn --deploy-mode cluster \
    --driver-memory 4g \
    --num-executors 4 --executor-cores 8 --executor-memory 30g \
    --conf spark.yarn.executor.memoryOverhead=600 
    --conf spark.local.dir=/home/hdfs/spark"


date
# Align the reads
#${dir_bin}/bwa mem -K 100000000 -v 3 -t $threads -R"$ReadGroup" $hg19  $R1 $R2 | ${dir_bin}/samtools view -@ $threads -S -b > $dir_analyzed_sample/$sp_name-bwa.bam
${dir_bin}/bwa mem -K 100000000 -v 3 -t $threads -R"$ReadGroup" $hg19  $R1 $R2 | ${dir_bin}/samtools view -@ 6 -1 -h  -F 4 -S -b > $dir_analyzed_sample/$sp_name-ali.bam
date

# filter reads un-mapped reads and sort bam file
#time bash -c "taskset -c 0-17 ${dir_bin}/samtools view -1 -@ $cores -b -h  -F 4 $dir_analyzed_sample/$sp_name-bwa.bam  > $dir_analyzed_sample/$sp_name-ali.bam"

# Sorting BAM
time bash -c "taskset -c 0-17  ${dir_bin}/samtools sort -@ $cores $dir_analyzed_sample/$sp_name-ali.bam  -o $dir_analyzed_sample/$sp_name-ali-sorted.bam"

# Don't need it, as we add the reading group in the bwa mem
# Add ReadGroup
#time bash -c "taskset -c 0-17  java -jar ${dir_bin}/picard.jar AddOrReplaceReadGroups  I=$dir_analyzed_sample/$sp_name-ali-sorted.bam  O=$dir_analyzed_sample/$sp_name-ali-sorted-RG.bam  RGSM=$sp_name RGLB=project RGPL=illumina RGID=none RGPU=none VALIDATION_STRINGENCY=LENIENT"


hdfs dfs -rm $sp_name-*
hdfs dfs -put $dir_analyzed_sample/$sp_name-ali-sorted.bam
date
# Remove duplicates
#time bash -c "gatk MarkDuplicatesSpark  --remove-all-duplicates false --create-output-bam-index --read-validation-stringency LENIENT -I hdfs://$hdfs_url/$sp_name-ali-sorted-RG.bam  -O hdfs://$hdfs_url/$sp_name-ali-sorted-RG-rmdup.bam  --metrics-file hdfs://$hdfs_url/$sp_name-pacard.metrics $spark_cmd"
#no metric for high perforamnce
gatk MarkDuplicatesSpark  --remove-all-duplicates false --create-output-bam-index --read-validation-stringency LENIENT -I hdfs://$hdfs_url/$sp_name-ali-sorted.bam  -O hdfs://$hdfs_url/$sp_name-ali-sorted-RG-rmdup.bam $spark_cmd


hdfs dfs -get $sp_name-ali-sorted-RG-rmdup.bam* $dir_analyzed_sample

date
# index BAM
#time bash -c "taskset -c 0-17  ${dir_bin}/samtools index -@ $cores $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam"

# RealignerTargetCreator
time bash -c "taskset -c 0-17  java -jar $dir_bin/GenomeAnalysisTK.jar -nt $cores -T RealignerTargetCreator -known $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf  -known $dir_vcf/1000G_phase1.indels.hg19.sites.vcf  -R $hg19 -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam -o $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam.list"

# IndelRealigner
time bash -c "taskset -c 0-17  java -jar $dir_bin/GenomeAnalysisTK.jar  -T IndelRealigner -known $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -known $dir_vcf/1000G_phase1.indels.hg19.sites.vcf -R $hg19 -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam -targetIntervals $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup.bam.list -o $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.bam"

#BaseRecalibrator
#time bash -c "taskset -c 0-17  java -jar $dir_bin/GenomeAnalysisTK.jar  -nct $cores -T BaseRecalibrator -l INFO  -R $hg19   -knownSites $dir_vcf/Mills_and_1000G_gold_standard.indels.hg19.sites.vcf -knownSites $dir_vcf/1000G_phase1.indels.hg19.sites.vcf  -I  $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.bam -o  $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.grp"
# GATK PrintReads
#time bash -c "taskset -c 0-17  java -jar $dir_bin/GenomeAnalysisTK.jar -nct 8 -T PrintReads -R $hg19 -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.bam -BQSR $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned.grp -o $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam"


#create script for Scatter Gather for bqsr and HaplotypeCaller
time python create_script.py
#BaseRecalibrator
time bash -c "taskset -c 0-17 parallel -j6 < ./script/bqsr.sh"
time ./script/merge_bqsr.sh
# GATK PrintReads
time bash -c "taskset -c 0-17 parallel -j6 < ./script/apply.sh"
time ./script/merge_bam.sh

# Filter mapping quality <10
#time bash -c "${dir_bin}/samtools view -@ $threads -h -q 10 $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal-0.bam |${dir_bin}/samtools view  -@ $threads -h -Sb > $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam"

# index BAM
time bash -c "taskset -c 0-17  ${dir_bin}/samtools index -@ $cores $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam"

# GATK UnifiedGenotyper
time bash -c "taskset -c 0-17  java -jar $dir_bin/GenomeAnalysisTK.jar -nt $cores -T UnifiedGenotyper -R $hg19  -L $bed  -metrics $dir_analyzed_sample/$sp_name-snps.metrics -stand_call_conf 10.0 -dcov 20000 -glm BOTH -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o $dir_analyzed_sample/$sp_name-Unified-SNP-INDLE.vcf"
#time java -jar $dir_bin/GenomeAnalysisTK.jar -nt $threads -T UnifiedGenotyper -R $hg19  -metrics $dir_analyzed_sample/$sp_name-snps.metrics -stand_call_conf 10.0 -stand_emit_conf 5.0 -dcov 20000 -glm BOTH -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o $dir_analyzed_sample/$sp_name-Unified-SNP-INDLE.vcf

#GARK HaplotypeCaller
#time java -jar $dir_bin/GenomeAnalysisTK.jar -nct  $threads -T HaplotypeCaller -R $hg19 -L $bed  -stand_call_conf 10.0 -stand_emit_conf 5.0  -minPruning 3  -mbq 5  -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o $dir_analyzed_sample/$sp_name-Haploy-SNP-INDLE.vcf
#time bash -c "taskset -c 0-17 java -jar $dir_bin/GenomeAnalysisTK.jar -nct 4 -T HaplotypeCaller -R $hg19 -stand_call_conf 10.0 -minPruning 3  -mbq 5  -I $dir_analyzed_sample/$sp_name-ali-sorted-RG-rmdup-realigned-recal.bam -o $dir_analyzed_sample/$sp_name-Haploy-SNP-INDLE.vcf"
time bash -c "taskset -c 0-17 parallel -j6 < ./script/htc.sh"
time ./script/merge_htc.sh

