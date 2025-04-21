#!/usr/bin/env bash

#SBATCH --job-name  filter_cells_rna
#SBATCH --partition LocalQ
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 13
#SBATCH --array=1-20
#SBATCH --open-mode=append
#SBATCH --output logs/peaks.log

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