#!/usr/bin/env Rscript --vanilla

####################################################################################################
####################################################################################################
#Filename:       filter_cells_rna.R
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


# IO paths
in.path <- file.path("processed_datasets",PROJ,"cellranger_output", paste0(ID,"-RNA"),"outs","filtered_feature_bc_matrix")
out.path <- file.path("processed_datasets",PROJ,"cell_labels","doublet_labels")
fig.path <- file.path("processed_datasets",PROJ,"figures")


# computeDoubletDensity method
# https://bioconductor.org/packages/3.15/bioc/vignettes/scDblFinder/inst/doc/computeDoubletDensity.html
library(DropletUtils)
sce <- read10xCounts(in.path)
sce

library(scuttle)
set.seed(1001)
sce <- logNormCounts(sce)

library(scran)
dec <- modelGeneVar(sce)
hvgs <- getTopHVGs(dec)

library(scater)
set.seed(1002)
sce <- runPCA(sce, ncomponents=30, subset_row=hvgs)
sce <- runTSNE(sce, dimred="PCA")

set.seed(1003)
library(scDblFinder)
scores <- computeDoubletDensity(sce, subset.row=hvgs)
cells <- colData(sce)
cells <- cells$Barcode
names(scores) <- cells
p1 <- plotTSNE(sce, colour_by=I(log1p(scores)))

# save figures and output
ggsave(filename = file.path(fig.path, paste0(ID,"_doublets_tsne_RNA.tiff")),plot = p1, width = 12,height = 8, device = "tiff")
ggsave(filename = file.path(fig.path, paste0(ID,"_doublets_tsne_RNA.jpg")),plot = p1, width = 12,height = 8, device = "jpg")
saveRDS(scores,file.path(out.path,paste0(ID,"_doublet_scores.rds")))

doublet_data <- readRDS(file = file.path(ann.path,paste0(ID,"_doublet_scores.rds")))
doublet_threshold <- 1.9


data <- readRDS(file.path(in.path,paste0(ID,"_Initial_Seurat_Object.rds")))
data <- AddMetaData(
  object = data,
  metadata = doublet_data,
  col.name = 'doublets'
)
filtered <- subset(
  x = data,
  subset = nFeature_RNA > 200 &
    percent.mt < 20 &
    doublets <= doublet_threshold
)

#######################
#  Output             #
#######################

saveRDS(filtered,file.path(out.path,paste0(ID,"_Filtered_Seurat_Object_RNA.rds")))

toc()
cat(Sys.info())
var$EXIT_STATUS <- 1
write(var$EXIT_STATUS, stdout())
