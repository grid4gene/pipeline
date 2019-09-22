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
export bg37="$base/reference/hs37d5.fa"
export dir_vcf="$base/reference"
export threads=26
export cores=26
export sp_name="deepvariant-bg37-$sample"
export dir_analyzed_sample="${sp_dir}/output"
export dir_bin='/home/test/WGS_pipeline/TOOLS/bin'


# Align the reads
time bash -c " ${dir_bin}/bwa mem  -t $threads  $bg37  $R1 $R2 | ${dir_bin}/samtools view -@ $threads -S -b > $dir_analyzed_sample/$sp_name-bwa.bam"
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