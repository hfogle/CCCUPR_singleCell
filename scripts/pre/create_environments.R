####################
# Conda Environments
####################

library(reticulate)

conda_create(
  envname = "singlecell",
  packages = c("mamba","pip","jq","nano","fastqc","seqfu","multiqc"),
  forge = TRUE,
  channel = c("bioconda","conda-forge","defaults"),
  conda = "/home/hfogle/miniconda3/bin/conda",
  python_version = "3.12")

conda_export(
  "singlecell",
  file = file.path(getwd(),"config","singlecell_conda_environment.yml"),
  json = FALSE,
  conda = "auto"
)