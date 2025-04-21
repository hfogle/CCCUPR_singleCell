#!/usr/bin/env Rscript --vanilla

####################################################################################################
####################################################################################################
#Filename:       label_cells_rna.R
#Author:         Homer Fogle, hfogle@cccupr.org
#Affiliation:    www.thefrancolab.org
#Last Updated:   21 Apr 2025

#Description:  
#
#
#
####################################################################################################
####################################################################################################
#######################
#  Required Packages  #
#######################
suppressPackageStartupMessages(library(argparse))   # Parse command line options
suppressPackageStartupMessages(library(tictoc))     # Runtime
suppressPackageStartupMessages(library(ggplot2))    # Plot generation
suppressPackageStartupMessages(library(future))     # Multicore settings
suppressPackageStartupMessages(library(jsonlite))   # Parse metadata files
suppressPackageStartupMessages(library(dplyr))      # Data manipulation
suppressPackageStartupMessages(library(tidyr))      # Data manipulation

source("scripts/functions.R")                     # Common script tasks

#######################
#  Script Parameters  #
#######################
var <- setvars(args = commandArgs(trailingOnly = TRUE))
dir <- setdirs(study = var$STUDY)
var$parser <- ArgumentParser()
var$parser$add_argument("-c", "--cores", type="integer", default=1, help="Parallel CPU cores per task [default %(default)s]",metavar="number")
var$parser$add_argument("-m", "--mem", type="integer", default=20, help="Memory in Gb per task [default %(default)s]",metavar="number")
var$parser$add_argument("-g", "--gpu", action="store_true", default=FALSE, help="Run on GPU node [default %(default)s]",metavar="number")
var$parser$add_argument("-i", "--id", type="character", default="", help="Parallel CPU cores per task [default %(default)s]",metavar="number")
var$parser$add_argument("-s", "--study", type="character", default="", help="Path to study working directory (e.g. ~/proj/CCCUPR_singleCell/data/study_id)) [default %(default)s]",metavar="number")
sink(file = dir$logs, append = FALSE, type = c("output", "message"))
var$TASK_MEM <- 20
var$TASK_CORES <- 4
var$EXIT_STATUS <- 0


#######################
#  Input              #
#######################
tic("RUNTIME: ")




#######################
#  Output             #
#######################
toc()
cat(Sys.info())
var$EXIT_STATUS <- 1
write(var$EXIT_STATUS, stdout())