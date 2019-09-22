#!/bin/bash
# Copyright 2018 Google LLC.

set -euo pipefail

## Preliminaries
# Set a number of shell variables, to make what follows easier to read.
BIN_VERSION="0.8.0"

BASE="/mnt/disk_nvme/20190826A1B1-1"
INPUT_DIR="${BASE}/input"
OUTPUT_DIR="${BASE}/output"

REF_DIR="/home/reference"
REF="hs37d5.fa.gz"

N_SHARDS="52"

#SP_NAME="S065_jiade-A_HGC-1911366_AHMYKWDSXX_S1_L002"
SP_NAME="S066_jiade-A_HGC-1911367_AHMYKWDSXX_S2_L002"

BAM="deepvariant-bg37-${SP_NAME}-ali-sorted-RG-rmdup.bam"

EXAMPLES="${OUTPUT_DIR}/${SP_NAME}.examples.tfrecord@${N_SHARDS}.gz"
GVCF_TFRECORDS="${OUTPUT_DIR}/${SP_NAME}.gvcf.tfrecord@${N_SHARDS}.gz"
CALL_VARIANTS_OUTPUT="${OUTPUT_DIR}/${SP_NAME}.cvo.tfrecord.gz"
OUTPUT_VCF="${SP_NAME}.output.vcf.gz"
OUTPUT_GVCF="${SP_NAME}.output.g.vcf.gz"
LOG_DIR="${OUTPUT_DIR}/logs"

CAPTURE_BED="agilent_sureselect_human_all_exon_v5_b37_targets.bed"

## Create local directory structure
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${INPUT_DIR}"
mkdir -p "${LOG_DIR}"


echo "Run DeepVariant..."
#  --regions="/reference/${CAPTURE_BED}" \
sudo docker run \
  -v "${REF_DIR}":"/reference" \
  -v "${INPUT_DIR}":"/input" \
  -v "${OUTPUT_DIR}:/output" \
  gcr.io/deepvariant-docker/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant \
  --model_type=WES \
  --ref="/reference/${REF}" \
  --reads="/output/${BAM}" \
  --output_vcf="/output/${OUTPUT_VCF}" \
  --output_gvcf="/output/${OUTPUT_GVCF}" \
  --num_shards=${N_SHARDS}
echo "Done."
echo


