# Opt-in native FACETS-CH integration test on synthetic counts and builtin hg38 GC.
# Does not validate capture/PoN suitability or establish mosaic sensitivity.
source("facets_grch38.R")
source("facets_segments.R")
stopifnot(requireNamespace("facetsCH",quietly=TRUE),requireNamespace("pctGCdata",quietly=TRUE))
n <- 600L
i <- seq_len(n)
normal_depth <- as.integer(120+i%%19)
loss <- i<=150L
sample_depth <- as.integer(round(normal_depth * ifelse(loss,0.75,1) * (1+(i%%7-3)/100)))
counts <- data.frame(Chromosome="chr1",Position=1000000L+i*1000L,
    NOR.DP=normal_depth,NOR.RD=as.integer(round(normal_depth/2)),
    TUM.DP=sample_depth,
    TUM.RD=as.integer(round(sample_depth*ifelse(loss,2/3,0.5))),
    catalog_kind="SNP",masked=FALSE)
bias <- data.frame(snp_index=paste("1",counts$Position,sep="_"),
    correction_bias=0,bias_unit="ln_alt_ref",allelic_pseudocount=1/6)
config <- list(reference_build="GRCh38",chromosome_map=c(chr1="1"),
    ndepth=20,ndepthmax=2000,het_threshold=0.1,snp_neighborhood=0,
    initial_cval=25,final_cval=50,gamma=0.25,min_nhet=15,seed=0,
    allelic_pseudocount=1/6,sample_het_min_depth=50,het_scale=TRUE,delta_cn=0)
search_before <- search()
result <- run_grch38_joint(counts,bias,config)
stopifnot(identical(search_before,search()))
stopifnot(result$status=="JOINT_SEGMENTATION_CANDIDATES",
    nrow(result$locus_audit)==n,
    all(is.finite(result$signal$cnlr)),
    all(is.finite(result$signal$valor)),
    all(result$locus_audit$allelic_status=="CORRECTED_ALLELIC_SIGNAL"),
    !result$sensitivity_calibrated,!result$copy_neutral_loh_confirmed)
signal <- result$fit$jointseg
stopifnot(nrow(signal)>0,length(unique(signal$seg))>=2)
segments <- summarize_native_segments(signal,config$chromosome_map)
stopifnot(sum(segments$n_markers)==nrow(signal),all(is.na(segments$evaluable_bp)))
cat("NATIVE_FACETS_GRCH38_SMOKE_PASS\n")
cat("facetsCH=",as.character(packageVersion("facetsCH")),
    "; pctGCdata=",as.character(packageVersion("pctGCdata")),"\n",sep="")
cat("native_engine_tested=TRUE; mosaic_sensitivity_measured=FALSE\n")
