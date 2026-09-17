#!/usr/bin/env Rscript
# Local native entrypoint. No upstream legacy CLI or implicit best-normal choice.
args <- commandArgs(trailingOnly=TRUE)
required <- c("--matrix","--bias","--config","--chromosome-map","--output","--library")
if (identical(args,"--help")) {
    cat("Usage: Rscript --vanilla facets_run.R",
        paste(paste(required,"PATH"),collapse=" "),"\n")
    quit(status=0)
}
if (length(args)!=2*length(required)) stop("exact named path arguments required; see --help")
keys <- args[seq(1,length(args),2)]
values <- args[seq(2,length(args),2)]
if (anyDuplicated(keys) || !setequal(keys,required) || any(!nzchar(values)))
    stop("invalid/duplicate command arguments")
options <- setNames(as.list(values),keys)
for (key in setdiff(required,"--output")) {
    if (!file.exists(options[[key]])) stop(paste("missing input",key))
    options[[key]] <- normalizePath(options[[key]],mustWork=TRUE)
}
library_path <- options[["--library"]]
if (!dir.exists(library_path) || !dir.exists(file.path(library_path,"facetsCH")))
    stop("facetsCH must be installed in the explicit isolated library")
# Suppress user/site libraries. All non-base dependencies must be installed
# in the declared library; the system R base/recommended library remains.
.libPaths(c(library_path,.Library),include.site=FALSE)
script_arg <- grep("^--file=",commandArgs(trailingOnly=FALSE),value=TRUE)
if (length(script_arg)!=1) stop("cannot resolve entrypoint location")
script_dir <- dirname(normalizePath(sub("^--file=","",script_arg),mustWork=TRUE))
source(file.path(script_dir,"facets_grch38.R"))
source(file.path(script_dir,"facets_segments.R"))
source(file.path(script_dir,"facets_io.R"))
run_joint_files(options[["--matrix"]],options[["--bias"]],options[["--config"]],
                options[["--chromosome-map"]],options[["--output"]])
