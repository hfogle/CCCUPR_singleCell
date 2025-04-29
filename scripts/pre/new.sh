#!/usr/bin/env bash
#SBATCH --job-name=new
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --output new_study_%j.log


# Activate Conda environment
source activate singlecell
echo "\n########################################################################"
############################
# Script Argument Parsing  #
#                          #
############################

Help()
{
   # Display Help
   echo "Add description of the script functions here."
   echo
   echo "Syntax: scriptTemplate [-s|o|h]"
   echo "options:"
   echo "s     Path to studydesign.json meta data file."
   echo "o     Parent directory for pipeline outputs. Must be writable"
   echo "h     Help options."
   echo
}

while getopts ":s:o:h" option; do
   case $option in
      h) # display Help
         Help
         exit;;
      s) # Enter metadata file path
         STUDYDESIGN=$(realpath $OPTARG)
         echo $STUDYDESIGN
         ;;
      o) # Enter output directory
         PROJ_DIR=$(realpath $OPTARG)
         echo $PROJ_DIR
         ;;
     \?) # Invalid option
         echo "Error: Invalid option"
         exit;;
   esac
done




############################
# Create study directories #
# Copy metadata files.     #
############################
echo "\n########################################################################"

if [! -w ${PROJ_DIR}/ ];
then
PROJ_DIR=${PWD}/data
echo "Defaulting to ${PWD}/data for study output."
fi

STUDY_ID=$(jq '.record.study_label' $STUDYDESIGN | tr -d \")
OUT_DIR=${PROJ_DIR}/${STUDY_ID}


if [ -w ${PROJ_DIR}/ ];
then
echo "Creating work directories for new study....."
mkdir $OUT_DIR
mkdir ${OUT_DIR}/reference_data
mkdir ${OUT_DIR}/reference_data/cellranger
mkdir ${OUT_DIR}/reference_data/inferCNV
mkdir ${OUT_DIR}/reference_data/scimilarity
mkdir ${OUT_DIR}/reference_data/cellcycle
mkdir ${OUT_DIR}/reference_data/ontologies
mkdir ${OUT_DIR}/reference_data/velocyto

mkdir ${OUT_DIR}/meta_data
mkdir ${OUT_DIR}/exchange_data
mkdir ${OUT_DIR}/raw_data
mkdir ${OUT_DIR}/reports
mkdir ${OUT_DIR}/figures
mkdir ${OUT_DIR}/logs

mkdir ${OUT_DIR}/processed_data
mkdir ${OUT_DIR}/processed_data/sample_alignment
mkdir ${OUT_DIR}/processed_data/sample_count
mkdir ${OUT_DIR}/processed_data/sample_filtered
mkdir ${OUT_DIR}/processed_data/sample_labeled
mkdir ${OUT_DIR}/processed_data/sample_velocity
mkdir ${OUT_DIR}/processed_data/cell_labels
mkdir ${OUT_DIR}/processed_data/integrated
mkdir ${OUT_DIR}/processed_data/sample_cnv
mkdir ${OUT_DIR}/analysis
mkdir ${OUT_DIR}/reports/fastqc
mkdir ${OUT_DIR}/reports/cellranger
mkdir ${OUT_DIR}/reports/checksums
mkdir ${OUT_DIR}/reports/seqfu
else
echo "NOT WRITABLE, Change Permissions"
fi

META_DIR=$(jq '.data.metadata_path' $STUDYDESIGN | tr -d \")
if [ ! -f ${META_DIR}/studydesign.json ]; then
   echo "Copying studydesign.json metadata file to working directory."
   cp ${META_DIR}/studydesign.json ${OUT_DIR}/meta_data/
else
   echo "File ${META_DIR}/studydesign.json does not exist."
fi
if [ ! -f ${META_DIR}/samplesheet.tsv ]; then
   echo "Copying samplesheet.tsv metadata file to working directory."
   cp ${META_DIR}/samplesheet.tsv ${OUT_DIR}/meta_data/
else
   echo "File ${META_DIR}/samplesheet.tsv does not exist."
fi

############################
# Raw data integrity checks#
# Perform checksums        #
############################
echo "\n########################################################################"

RAW_DIR=$(jq '.data.raw_path' $STUDYDESIGN | tr -d \")


### create symlinks to raw files in working directory
echo "Creating raw data symlinks in working directory....."
for file in $(ls ${RAW_DIR}/*fastq.gz)
do
base=$(basename "$file")
ln -s ${file} ${OUT_DIR}/raw_data/${base}
done

### sha256 checksums
for file in $(ls ${RAW_DIR}/*fastq.gz) 
do
if [ -f ${file}.sha256 ] 
then
sum=$(sha256sum --check ${file}.sha256 $file)
echo "Checksum status of $file ..... $sum"
else
echo "Creating sha256 checksum for $file"
base=$(basename "$file")
sha256sum $file > ${OUT_DIR}/reports/checksums/${base}.sha256
fi
done
echo "\n########################################################################"

### Count raw read files
R1_COUNT=$(ls -1q ${RAW_DIR}/*R1_001.fastq.gz | wc -l)
R2_COUNT=$(ls -1q ${RAW_DIR}/*R2_001.fastq.gz | wc -l)
echo "R1 File Count" $R1_COUNT
echo "R2 File Count" $R2_COUNT
echo "\n########################################################################"

### Validate fastq files (SeqFu)
echo "Running SeqFu FASTQ file integrity checks ....."
seqfu check --thousands --dir $RAW_DIR > ${OUT_DIR}/reports/seqfu/fastq_file_integrity_report.tsv
echo "\n########################################################################"

### Generate FastQC reports
for file in $(ls ${RAW_DIR}/*fastq.gz) 
do
fastqc --noextract --nogroup -o ${OUT_DIR}/reports/fastqc $file
done
multiqc -o ${OUT_DIR}/reports -n fastqc_summary_report.html ${OUT_DIR}/reports/fastqc 
############################
# Redirect Log             #
#                          #
############################

mv new_study_$SLURM_JOB_ID.log ${OUT_DIR}/logs/