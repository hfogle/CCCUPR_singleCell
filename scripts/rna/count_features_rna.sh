#!/usr/bin/env bash
#SBATCH --job-name  count_features_rna
#SBATCH --partition LocalQ
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 17
#SBATCH --array=1-$2


WDIR=$pwd
STUDY=$1
STUDYDESIGN=${WDIR}/data/${STUDY}/meta_data/studydesign.json
SAMPLESHEET=${WDIR}/data/${STUDY}/meta_data/samplesheet.tsv

ATICOL=awk -v RS='\t' '/^ARRAY_ID$/{print NR; exit}' $SAMPLESHEET
SRCCOL=awk -v RS='\t' '/^SOURCE_ID$/{print NR; exit}' $SAMPLESHEET
BATCOL=awk -v RS='\t' '/^BATCH_ID$/{print NR; exit}' $SAMPLESHEET
EXTCOL=awk -v RS='\t' '/^EXTRACT_ID$/{print NR; exit}' $SAMPLESHEET
FIDCOL=awk -v RS='\t' '/^FILE_ID$/{print NR; exit}' $SAMPLESHEET
SAMCOL=awk -v RS='\t' '/^SAMPLE_ID$/{print NR; exit}' $SAMPLESHEET
GRPCOL=awk -v RS='\t' '/^GROUP_ID$/{print NR; exit}' $SAMPLESHEET
SAMPLE=$(awk -v ARRAY_ID=$SLURM_ARRAY_TASK_ID '$ATICOL==ARRAY_ID {print $SAMCOL}' $SAMPLESHEET)
FILEID=$(awk -v ARRAY_ID=$SLURM_ARRAY_TASK_ID '$ATICOL==ARRAY_ID {print $FIDCOL}' $SAMPLESHEET)
GROUP=$(awk -v ARRAY_ID=$SLURM_ARRAY_TASK_ID '$ATICOL==ARRAY_ID {print $GRPCOL}' $SAMPLESHEET)

logs=${WDIR}/data/${STUDY}/processed_data/logs
outdir=${WDIR}/data/${STUDY}


REF=$(jq '.data.alignment_reference_path' $STUDYDESIGN | tr -d \")
RAW_DIR=$(jq '.data.raw_path' $STUDYDESIGN | tr -d \")

cellranger count --id=${SAMPLE} \
           --transcriptome=${REF} \
           --fastqs=${RAW_DIR} \
           --sample=${SAMPLE} \
           --localcores=16 \
           --localmem=64 \
           --create-bam=true \
           --output-dir=${outdir}/processed_data/sample_alignment/${SAMPLE}-RNA &> ${logs}/cellranger_count_RNA_${SAMPLE}.log


cp ${outdir}/processed_data/sample_alignment/${SAMPLE}-RNA/outs/web_summary.html ${outdir}/reports/cellranger/${SAMPLE}_cellranger_report.html

echo "########################################################"
