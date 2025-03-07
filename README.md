# CCCUPR_singleCell
Data Processing Pipeline for scRNA-Seq, scATAC-Seq, Paired scRNA-scATAC-Seq, and Multiome scRNA-scATAC-Seq Datasets.

# Contents
- [Summary](#id-section1)
- [Installation](#id-section2)
- [Input Data Formatting](#id-section3)
- [Pre-Processing](#id-section4)
- [scRNA-Seq Only Workflow](#id-section5)
- [scATAC-Seq Only Workflow](#id-section6)
- [Paired scRNA-scATAC-Seq Workflow](#id-section7)
- [Multiome scRNA-scATAC-Seq Workflow](#id-section8)
- [Output Data Structure](#id-section9)
- [Documentation](#id-section10)

<div id='id-section1'/>
## Summary
<div id='id-section2'/>
## Installation
<div id='id-section3'/>
## Input Data Formatting
<div id='id-section4'/>
## Pre-Processing
<div id='id-section5'/>
## scRNA-Seq Only Workflow
<div id='id-section6'/>
## scATAC-Seq Only Workflow
<div id='id-section7'/>
## Paired scRNA-scATAC-Seq Workflow
<div id='id-section8'/>
## Multiome scRNA-scATAC-Seq Workflow
<div id='id-section9'/>
## Output Data Structure
<div id='id-section10'/>
## Documentation

### Directory Structure

CCCUPR_singleCell
  - data
    - <STUDY_ID>
      - preprocessed_data
      - reference_data
      - processed_data
        - count
        - inital
        - filtered
        - labeled
          - cell_labels
        - integrated
      - meta_data
      - reports
      - figures
      - logs
  - scripts
    - pre
      - new.sh
      - get_genome_reference_sources.sh
      - make_custom_genome_reference.sh
      - build_genome_reference_data.sh
      - get_labeling_reference_sources.sh
      - build_labeling_references.sh
      - standardize_input_data.sh
      - install_libraries.R
    - rna
      - count_features_rna.sh
      - create_assay_rna.R
      - filter_cells_rna.sh
      - make_filter_labels_rna.R
      - label_cells_rna.sh
      - make_celltype_labels_rna.R
      - make_CNV_labels_rna.R
      - integrate_samples_rna.sh
      - integrate_samples_rna.R
    - atac
      - count_feature_atac.sh
      - create_assay_atac.R
      - filter_cells_atac.sh
      - label_cells_atac.sh
      - create_assay_activity.sh
      - integrate_samples_activity.sh
      - integrate_samples_atac.sh   
    - dual
      - feature_count_multi.sh
      - initialize_objects_multi.sh
      - filter_cells_multi.sh
      - label_cells_multi.sh
      - transfer_labels_paired.sh
      - integrate_samples_paired.sh 
    - analysis  
  - config
    - conda_environment.yml
    - CCCUPR_singleCell.dockerfile
    - meta_data_model.json
    - workflow_template.json
    - workflow.snakemake
  - docs
    - example_reports
    - samplesheet_template.tsv
    - studydesign_template.json
    - workflow_template.json
    - workflow_documentation.md
