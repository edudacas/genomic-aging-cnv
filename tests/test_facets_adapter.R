# Contract tests with an injected fake engine; NOT a native FACETS-CH test.
source("facets_grch38.R")
config <- list(reference_build="GRCh38", chromosome_map=c(chr1="1",chrY="Y"),
    ndepth=20, ndepthmax=2000, het_threshold=0.1, snp_neighborhood=0,
    initial_cval=25, final_cval=150, gamma=0.25, min_nhet=15, seed=0,
    allelic_pseudocount=1/6, sample_het_min_depth=50, het_scale=TRUE,delta_cn=0)
counts <- data.frame(Chromosome=c(rep("chr1",4),"chrY"),Position=c(20,30,40,50,70),
    NOR.DP=100,NOR.RD=50,TUM.DP=100,TUM.RD=50,
    catalog_kind=c("SNP","DEPTH","SNP","SNP","SNP"),masked=c(FALSE,FALSE,FALSE,TRUE,FALSE))
counts$catalog_row <- c(2,5,10,20,30)
bias <- data.frame(snp_index="1_20",correction_bias=0.2,
                   bias_unit="ln_alt_ref",allelic_pseudocount=1/6)
engine <- new.env()
engine$procSnps <- function(x,ndepth,het.thresh,snp.nbhd,nX,unmatched,ndepthmax) {
    stopifnot(nX==23,unmatched,ndepth==20,ndepthmax==2000)
    data.frame(chrom=as.integer(x$Chromosome),maploc=x$Position,rCountT=x$TUM.DP)
}
engine$counts2logROR <- function(x,gbuild,unmatched) {
    stopifnot(gbuild=="hg38",unmatched)
    transform(x,cnlr=c(0,0.02,-0.01),valor=log(2),het=1)
}
engine$segsnps <- function(x,cval,hetscale,delta,gamma) {
    stopifnot(cval==25,hetscale,delta==0,gamma==0.25)
    stopifnot(abs(x$valor[1]-(log(2)-0.2))<1e-12)
    stopifnot(is.na(x$valor[2]),is.na(x$valor[3]),all(x$het[2:3]==0))
    list(seg.tree=list(),jointseg=x)
}
engine$procSample <- function(x,cval,min.nhet) {
    stopifnot(x$gbuild=="hg38",x$nX==23,cval==150,min.nhet==15)
    list(jointseg=x$jointseg)
}
result <- run_grch38_joint(counts,bias,config,engine)
stopifnot(result$status=="JOINT_SEGMENTATION_CANDIDATES",
    !result$copy_neutral_loh_confirmed,!result$somatic_origin_classified,
    !result$sensitivity_calibrated,nrow(result$exclusions)==2,
    "Y_NOT_SUPPORTED_BY_NATIVE_JOINT_ENGINE" %in% result$exclusions$reason,
    "FROZEN_RESOURCE_MASK" %in% result$exclusions$reason)
stopifnot(nrow(result$locus_audit)==nrow(counts),
    identical(result$locus_audit$allelic_status,
              c("CORRECTED_ALLELIC_SIGNAL","DEPTH_ONLY","BIAS_UNAVAILABLE",
                "NOT_EVALUATED","NOT_EVALUATED")),
    identical(result$signal$source_input_row,1:3))
stopifnot(identical(result$signal$source_catalog_row,c(2,5,10)),
          identical(result$locus_audit$catalog_row,counts$catalog_row))
must_fail <- function(expression) {
    result <- tryCatch({force(expression); FALSE},error=function(e) TRUE)
    stopifnot(result)
}
bad <- config; bad$reference_build <- "GRCh37"
must_fail(run_grch38_joint(counts,bias,bad,engine))
bad <- config; bad$sample_het_min_depth <- 20
must_fail(run_grch38_joint(counts,bias,bad,engine))
bad <- config; bad$allelic_pseudocount <- 0
must_fail(run_grch38_joint(counts,bias,bad,engine))
bad <- bias; bad$bias_unit <- "log2_alt_ref"
must_fail(run_grch38_joint(counts,bad,config,engine))
bad <- counts; bad$TUM.RD[1] <- 101
must_fail(run_grch38_joint(bad,bias,config,engine))
bad <- counts; bad$Chromosome[1] <- "1"
must_fail(run_grch38_joint(bad,bias,config,engine))
must_fail(run_grch38_joint(counts[c(2,1,3,4,5),],bias,config,engine))
bad <- bias; bad$correction_bias <- NA_real_
must_fail(run_grch38_joint(counts,bad,config,engine))
engine$counts2logROR <- function(x,gbuild,unmatched) transform(x,cnlr=0,valor=log(2),het=1)
must_fail(run_grch38_joint(counts,bias,config,engine))
cat("FACETS_GRCH38_ADAPTER_CONTRACT_PASS; native_engine_tested=FALSE\n")
