#!/usr/bin/env bash

#SBATCH --job-name  count_features_rna
#SBATCH --partition LocalQ
#SBATCH --ntasks 1
#SBATCH --cpus-per-task 13
#SBATCH --array=1-2
#SBATCH --open-mode=append
#SBATCH --output logs/count_features_rna.log

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


table=/opt/data/sc-RNA-seq/PBMC_demo_scRNA/meta_data/samplesheet.tsv
indir=/opt/data/sc-RNA-seq/PBMC_demo_scRNA/raw_data
outdir=/opt/data/sc-RNA-seq/PBMC_demo_scRNA/processed_data/sample_alignment
ref=/opt/data/sc-RNA-seq/PBMC_demo_scRNA/reference_data/refdata-gex-GRCh38-2024-A
logs=/opt/data/sc-RNA-seq/PBMC_demo_scRNA/logs

# Get sample ID from the samplesheet

sample=$(awk -v ARRAY_ID=$SLURM_ARRAY_TASK_ID '$1==ARRAY_ID {print $8}' $table)
prefix=$(awk -v ARRAY_ID=$SLURM_ARRAY_TASK_ID '$1==ARRAY_ID {print $7}' $table)



echo "########################################################"


cellranger count --id=${sample} \
           --transcriptome=${ref} \
           --fastqs=${indir} \
           --sample=${sample} \
           --localcores=16 \
           --localmem=64 \
           --create-bam=true \
           --output-dir=${outdir}/${sample}-RNA &> ${logs}/cellranger_count_RNA_${sample}.log



echo "########################################################"
