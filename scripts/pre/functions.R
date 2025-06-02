### R Functions

#source()

# return standard script paths
setdirs <- function(){
  dir <- list()
  dir$study <- "elsb"
  dir$wdir <- file.path(getwd(),"data",dir$study)
  dir$config <- file.path(getwd(),"config")
  dir$reference <- file.path(dir$wdir,"reference_data")
  dir$metadata <- file.path(dir$wdir,"meta_data")
  dir$logs <- file.path(dir$wdir,"logs") 
  dir$figures <- file.path(dir$wdir,"figures") 
  dir$reports <- file.path(dir$wdir,"reports") 
  dir$inpath <- file.path(dir$wdir,"processed_data") 
  dir$infile <- ""
  dir$outpath <- file.path(dir$wdir,"processed_data") 
  dir$outfile <- "" 
  return(dir)
}
#return standard script parameters
setvars <- function(x){
  var <- list()
  var$ARGS <- (commandArgs(TRUE))
  if(missing(var$ARGS[1])) {
    var$ID <- "44F0AL"
  } else {
    var$ID <- var$ARGS[1]
  }
  if(missing(var$ARGS[2])) {
    var$STUDY <- "elsb"
  } else {
    var$STUDY <- var$ARGS[2]
  }
  if(missing(var$ARGS[3])) {
    var$TASK_CORES <- 1
  } else {
    var$TASK_CORES <- var$ARGS[3]
  }
  if(missing(var$ARGS[4])) {
    var$TASK_MEM <- 20*1024^3 # 20GB
  } else {
    var$TASK_MEM <- as.numeric(var$ARGS[4]) * 1024^3
  }
  var$SYS_MEM <- system('free',intern=TRUE)[2]
  var$SYS_MEM <- as.numeric(strsplit(var$SYS_MEM, "\\s+")[[1]][7]) * 1000
  var$SYS_CORES <- parallel::detectCores(all.tests = FALSE, logical = TRUE)
  var$PROCESSES <- 1
  if(interactive()){
    options(future.globals.maxSize= 0.9*var$SYS_MEM)
    plan("sequential")
  }else{
    options(future.globals.maxSize= var$TASK_MEM)
    plan("multicore", workers = var$PROCESSES)  
  }
  
  set.seed(123)  
  return(var)
}
# run rna clustering sequence
# clust_rna
# # run atac clustering sequence
# clust_atac
# # run activity clustering sequence
# clust_act
# generate standard clustering figures
umap_figures <- function(seurat, dir, ID){
  plots <- list()
  dir.create(file.path(dir, ID), showWarnings = FALSE)
  # return a list object of plot objects and saves id.rds, id.jpg, id.tif to figures subdirectory
  # just list all metadata columns and make a plot for each
  # gzip the jpg and the tiff files
  # rds for future touch ups
}

# build standard boilerplate content for pipeline script files
boilerplate <- function(name, type, mem, cores, gpu = FALSE, R = FALSE, slurm = FALSE, nextflow = FALSE){
  
  if (R == TRUE){
    file=file.path(getwd(),"scripts",paste0(name,".R"))
    cat("Rscript file created: ",file)
    sink(file = file, append = FALSE, type = "output")
    #######################
    # FILE CONTENT STARTS #
    #######################
    cat("#!/usr/bin/env Rscript --vanilla\n")
    cat("\n####################################################################################################")
    cat("\n####################################################################################################")
    cat("\n#Filename:      ",paste0(name,".R"))
    cat("\n#Author:        ", "Homer Fogle, hfogle@cccupr.org")
    cat("\n#Affiliation:   ", "www.thefrancolab.org")
    cat("\n#Last Updated:  ", format(Sys.Date(), "%d %b %Y"))
    cat("\n\n#Description: ", "\n#\n#\n#")
    cat("\n####################################################################################################")
    cat("\n####################################################################################################")
    
    cat("\n#######################")
    cat("\n#  Required Packages  #")
    cat("\n#######################")
    cat("\nsuppressPackageStartupMessages(library(argparse))   # Parse command line options")
    cat("\nsuppressPackageStartupMessages(library(tictoc))     # Runtime")
    cat("\nsuppressPackageStartupMessages(library(ggplot2))    # Plot generation")
    cat("\nsuppressPackageStartupMessages(library(future))     # Multicore settings")
    cat("\nsuppressPackageStartupMessages(library(jsonlite))   # Parse metadata files")
    cat("\nsuppressPackageStartupMessages(library(dplyr))      # Data manipulation")
    cat("\nsuppressPackageStartupMessages(library(tidyr))      # Data manipulation\n")
    cat("\nsource(\"scripts/functions.R\")                     # Common script tasks\n" ) 
    
    cat("\n#######################")
    cat("\n#  Script Parameters  #")
    cat("\n#######################")
    cat("\nvar <- setvars(args = commandArgs(trailingOnly = TRUE))")
    cat("\ndir <- setdirs(study = var$STUDY)")
    cat("\nvar$parser <- ArgumentParser()")
    cat("\nvar$parser$add_argument(\"-c\", \"--cores\", type=\"integer\", default=1, help=\"Parallel CPU cores per task [default %(default)s]\",metavar=\"number\")")
    cat("\nvar$parser$add_argument(\"-m\", \"--mem\", type=\"integer\", default=20, help=\"Memory in Gb per task [default %(default)s]\",metavar=\"number\")")
    
    cat("\nvar$parser$add_argument(\"-g\", \"--gpu\", action=\"store_true\", default=FALSE, help=\"Run on GPU node [default %(default)s]\",metavar=\"number\")")
    cat("\nvar$parser$add_argument(\"-i\", \"--id\", type=\"character\", default=\"\", help=\"Parallel CPU cores per task [default %(default)s]\",metavar=\"number\")")
    cat("\nvar$parser$add_argument(\"-s\", \"--study\", type=\"character\", default=\"\", help=\"Path to study working directory (e.g. ~/proj/CCCUPR_singleCell/data/study_id)) [default %(default)s]\",metavar=\"number\")")
    
    
    cat("\nsink(file = dir$logs, append = FALSE, type = c(\"output\", \"message\"))")
    cat("\nvar$TASK_MEM <- 20")
    cat("\nvar$TASK_CORES <- 4")
    cat("\nvar$EXIT_STATUS <- 0")
    cat("\n\n\n#######################")
    cat("\n#  Input              #")
    cat("\n#######################")
    cat("\ntic(\"RUNTIME: \")")
    cat("\n\n\n\n\n#######################")
    cat("\n#  Output             #")
    cat("\n#######################")
    
    #######################
    # FILE CONTENT ENDS   #
    #######################
    
    cat("\ntoc()")
    cat("\ncat(Sys.info())")
    cat("\nvar$EXIT_STATUS <- 1")
    cat("\nwrite(var$EXIT_STATUS, stdout())")
    sink()
    
  }
  
  if (slurm == TRUE){
    file=file.path("scripts",paste0(name,".sh"))
    sink(file = file, append = FALSE, type = "output")
    #######################
    # FILE CONTENT STARTS #
    #######################
    cat("#!/usr/bin/env bash\n")
    cat("\n#SBATCH --job-name ",name)
    cat("\n#SBATCH --partition LocalQ")
    cat("\n#SBATCH --ntasks 1")
    cat("\n#SBATCH --cpus-per-task 13")
    cat("\n#SBATCH --array=1-20")
    cat("\n#SBATCH --open-mode=append")
    cat("\n#SBATCH --output logs/peaks.log")
    cat("\n")
    
    
    cat("\nWDIR=$pwd")
    cat("\nSTUDY=$1")
    cat("\nSTUDYDESIGN=${WDIR}/data/${STUDY}/meta_data/studydesign.json")
    cat("\nSAMPLESHEET=${WDIR}/data/${STUDY}/meta_data/samplesheet.tsv")
    
    cat("\nATICOL=awk -v RS='\\t' '/^ARRAY_ID$/{print NR; exit}' $SAMPLESHEET")
    cat("\nSRCCOL=awk -v RS='\\t' '/^SOURCE_ID$/{print NR; exit}' $SAMPLESHEET")
    cat("\nBATCOL=awk -v RS='\\t' '/^BATCH_ID$/{print NR; exit}' $SAMPLESHEET")
    cat("\nEXTCOL=awk -v RS='\\t' '/^EXTRACT_ID$/{print NR; exit}' $SAMPLESHEET")
    cat("\nFIDCOL=awk -v RS='\\t' '/^FILE_ID$/{print NR; exit}' $SAMPLESHEET")
    cat("\nSAMCOL=awk -v RS='\\t' '/^SAMPLE_ID$/{print NR; exit}' $SAMPLESHEET")
    cat("\nGRPCOL=awk -v RS='\\t' '/^GROUP_ID$/{print NR; exit}' $SAMPLESHEET")
    
    cat("\nSAMPLE=$(awk -v ARRAY_ID=$SLURM_ARRAY_TASK_ID '$ATICOL==ARRAY_ID {print $SAMCOL}' $SAMPLESHEET)")
    cat("\nFILEID=$(awk -v ARRAY_ID=$SLURM_ARRAY_TASK_ID '$ATICOL==ARRAY_ID {print $FIDCOL}' $SAMPLESHEET)")
    cat("\nGROUP=$(awk -v ARRAY_ID=$SLURM_ARRAY_TASK_ID '$ATICOL==ARRAY_ID {print $GRPCOL}' $SAMPLESHEET)")
    cat("\n\nlogs=${WDIR}/data/${STUDY}/processed_data/logs")
    #######################
    # FILE CONTENT ENDS   #
    #######################
    sink()
  }
  
  if (nextflow == TRUE){
    file=file.path(get(),"scripts",paste0(name,".nextflow.txt"))
    sink(file = file, append = TRUE, type = "output")
    #######################
    # FILE CONTENT STARTS #
    #######################
    
    
    #######################
    # FILE CONTENT ENDS   #
    #######################
    sink()
  }
  
}