####################
# Install Packages
####################


install.packages("argparse")
install.packages("tidyverse")
install.packages("tictoc")
install.packages("future")
install.packages("RSQLite")
install.packages("reticulate")
install.packages("BiocManager")
install.packages("remotes")
install.packages("devtools")
install.packages("ggvenn")
install.packages("R.utils")

### Seurat Dependencies

remotes::install_github("satijalab/seurat", "seurat5", quiet = TRUE)
remotes::install_github("satijalab/seurat-data", "seurat5", quiet = TRUE)
devtools::install_github('cole-trapnell-lab/leidenbase')
remotes::install_github("satijalab/seurat-wrappers", "seurat5")
remotes::install_github("stuart-lab/signac", "seurat5", quiet = TRUE)
remotes::install_github("bnprks/BPCells/r")

### scCustomize Dependencies

install.packages("ggpubr")
install.packages("hdf5r")
devtools::install_version("rliger",version="2.1-0")
BiocManager::install("ComplexHeatmap")
BiocManager::install("dittoSeq")
BiocManager::install("DropletUtils")
BiocManager::install("Nebulosa")
install.packages("scCustomize")

### SuperCell Dependencies
install.packages("igraph")
install.packages("RANN")
install.packages("WeightedCluster")
install.packages("corpcor")
install.packages("weights")
install.packages("Hmisc")
install.packages("Matrix")
install.packages("patchwork")
install.packages("plyr")
install.packages("irlba")
remotes::install_github("GfellerLab/SuperCell") ### !!! github issue

### Escape Dependencies


BiocManager::install("BSgenome")
BiocManager::install("Biostrings")
BiocManager::install("org.Hs.eg.db")

BiocManager::install(c('BSgenome.Hsapiens.UCSC.hg38', 'EnsDb.Hsapiens.v86'))
BiocManager::install("EnhancedVolcano")
# BiocManager::install("genefu")
renv::install("genefu_2.40.0.tar.gz")


### scTyper
install.packages("fastqcr")
install.packages("parallel")
install.packages("kableExtra")
install.packages("gProfileR")
install.packages("perm")
install.packages("reshape2")
BiocManager::install("infercnv")
devtools::install_github("omicsCore/scTyper")
