args <- commandArgs(trailingOnly=TRUE)
stopifnot(length(args)==2)
library_dir <- normalizePath(args[1], mustWork=TRUE)
archives <- normalizePath(args[2], mustWork=TRUE)
.libPaths(c(library_dir, .Library))
options(repos=c(CRAN='https://cloud.r-project.org'), timeout=600)
required <- c('dplyr','data.table','ggplot2','argparse','stringr','patchwork','glue','R.utils')
available <- available.packages()
dependencies <- unique(c(required, unlist(tools::package_dependencies(
  required, available, which=c('Depends','Imports','LinkingTo'), recursive=TRUE))))
installed <- rownames(installed.packages())
missing <- setdiff(dependencies, installed)
if (length(missing)) {
  downloaded <- download.packages(missing, destdir=archives, available=available, type='source')
  # install.packages(repository) orders dependencies; retained source archives lock this attempt.
  install.packages(missing, lib=library_dir, type='source', Ncpus=2, destdir=archives)
}
stopifnot(all(vapply(required, requireNamespace, logical(1), quietly=TRUE)))
write.table(installed.packages(lib.loc=c(library_dir, .Library))[,c('Package','Version')],
            file=file.path(archives,'R-packages.tsv'), sep='\t', row.names=FALSE, quote=FALSE)
