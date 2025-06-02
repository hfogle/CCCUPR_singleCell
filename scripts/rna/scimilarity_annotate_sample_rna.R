#!/usr/bin/env Rscript --vanilla

####################################################################################################
####################################################################################################
#Filename:       scimilarity_annotate_sample_rna.R
#Author:         Homer Fogle, hfogle@cccupr.org
#Affiliation:    www.thefrancolab.org
#Last Updated:   01 May 2025

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

source("scripts/pre/functions.R")                     # Common script tasks

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

#!/usr/bin/env Rscript --vanilla

#### R Libraries
library(Seurat)
library(SeuratObject)
library(scCustomize)

#### Set Parameters
args=(commandArgs(TRUE))

var$blood.celltypes <- c("CD8-positive, alpha-beta T cell",
                         "CD4-positive, alpha-beta T cell",
                         "regulatory T cell",
                         "T cell",
                         "B cell",
                         "natural killer cell",
                         "dendritic cell",
                         "monocyte",
                         "platelet",
                         "myeloid cell",
                         "macrophage",
                         "hematopoietic stem cell",
                         "mast cell",
                         "eosiophil",
                         "neutrophil",
                         "basophil",
                         "erythrocyte")


var$breast.celltypes <- c("T cell",
                          "B cell",
                          "monocyte",
                          "macrophage",
                          "myeloid cell", 
                          "natural killer cell", 
                          "dentritic cell",
                          "epithelial cell",
                          "endothelial cell",
                          "basal cell",
                          "PVL",
                          "CAF",
                          "plasmablast")

dir$inpath <- file.path(dir$inpath,"sample_filtered")
dir$outpath <- file.path(dir$outpath, "cell_labels")
dir$reffile <- file.path(dir$reference,"scimilarity","model_v1.1")
dir$infile.seurat <- paste0(var$ID,"_Filtered_Seurat_Object_RNA.rds")
dir$infile.anndata <- paste0(var$ID,"_Filtered_AnnData_Object_RNA.h5ad")


#### Convert filtered data to anndata format for scimilarity
if (!file.exists(file.path(dir$inpath,dir$infile.anndata))){
  filtered.rna <- readRDS(file.path(dir$inpath,dir$infile.seurat))
  as.anndata(
    filtered.rna,
    file_path = dir$inpath,
    file_name = dir$infile.anndata,
    assay = "RNA",
    main_layer = "counts",
    other_layers = "data",
    transer_dimreduc = TRUE,
    verbose = TRUE
  )
  rm(filtered.rna)
}

#### Import dataset and celltype model
reticulate::use_condaenv("scimilarity",required = TRUE)
sc <- reticulate::import("scanpy")
data = sc$read(file.path(dir$inpath,dir$infile.anndata))
data$layers["counts"] = data$X

scimilarity <- reticulate::import("scimilarity")
ca = scimilarity$CellAnnotation(model_path=dir$refpath)

data = scimilarity$align_dataset(data, ca$gene_order)

data = scimilarity$lognorm_counts(data)

data$obsm["X_scimilarity"] = ca$get_embeddings(data$X)


pred = ca$get_predictions_knn(data$obsm["X_scimilarity"])
names(pred) <- c("predictions","nn_idxs","nn_dists","nn_stats")

data$obs["predictions_unconstrained"] = pred$predictions





obs <- data$obs
labels <- data.frame(predicted.id=as.character(obs$predictions_unconstrained))
rownames(labels) <- rownames(obs)



ca$safelist_celltypes(var$blood.celltypes)

data = ca$annotate_dataset(data)
labels["predictions_major"] <- data$obs["celltype_hint"]

saveRDS(labels, file = file.path(dir$outpath,paste0(var$ID,"_scimilarity_celltype_labels.rds")))


#######################
#  Output             #
#######################
library(DBI)
mydb <- dbConnect(RSQLite::SQLite(), "my-db.sqlite")
dbWriteTable(mydb, "mtcars", mtcars)
dbDisconnect(mydb)
toc()
cat(Sys.info())
var$EXIT_STATUS <- 1
write(var$EXIT_STATUS, stdout())