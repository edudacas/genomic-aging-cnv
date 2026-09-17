# GRCh38 adapter composing audited facetsCH primitives, not the upstream CLI.
# Native package pinned upstream in UPSTREAM_AUDIT.md. Caller installation,
# GC/resource provenance and input role matching are separate prerequisites.

require_fields <- function(x, fields, label) {
    if (!all(fields %in% names(x))) stop(paste(label, "missing required fields"))
}

run_grch38_joint <- function(rcmat, bias, config,
                             engine=asNamespace("facetsCH")) {
    require_fields(config, c("reference_build", "chromosome_map", "ndepth",
        "ndepthmax", "het_threshold", "snp_neighborhood", "initial_cval",
        "final_cval", "gamma", "min_nhet", "seed", "allelic_pseudocount",
        "sample_het_min_depth", "het_scale", "delta_cn"), "config")
    if (!identical(config$reference_build, "GRCh38"))
        stop("explicit GRCh38 required")
    numeric_fields <- c("ndepth", "ndepthmax", "het_threshold", "snp_neighborhood",
        "initial_cval", "final_cval", "gamma", "min_nhet", "seed",
        "allelic_pseudocount", "sample_het_min_depth", "delta_cn")
    for (field in numeric_fields) {
        value <- config[[field]]
        if (!is.numeric(value) || length(value)!=1L || !is.finite(value))
            stop(paste("invalid numeric config", field))
    }
    if (config$ndepth <= 0 || config$ndepthmax <= config$ndepth ||
        config$het_threshold <= 0 || config$het_threshold >= 0.5 ||
        config$snp_neighborhood < 0 || config$initial_cval <= 0 ||
        config$final_cval < config$initial_cval || config$gamma <= 0 ||
        config$min_nhet < 1 || config$delta_cn < 0 ||
        config$seed < 0 || config$seed > .Machine$integer.max ||
        config$seed != floor(config$seed))
        stop("invalid segmentation parameter range")
    if (!is.logical(config$het_scale) || length(config$het_scale)!=1L ||
        is.na(config$het_scale)) stop("explicit het_scale required")
    # These are hardcoded in the audited primitive. Never silently ignore a
    # requested alternative; a different native implementation needs new tests.
    if (config$allelic_pseudocount != 1/6 || config$sample_het_min_depth != 50)
        stop("requested estimator differs from audited native primitive")
    mapping <- config$chromosome_map
    if (!is.character(mapping) || is.null(names(mapping)) ||
        anyNA(mapping) || anyNA(names(mapping)) || any(names(mapping)=="") ||
        anyDuplicated(names(mapping)) || anyDuplicated(unname(mapping)) ||
        !all(mapping %in% c(as.character(1:22), "X", "Y")))
        stop("explicit one-to-one chromosome map required")
    fields <- c("Chromosome", "Position", "NOR.DP", "NOR.RD", "TUM.DP", "TUM.RD",
                "catalog_kind", "masked")
    require_fields(rcmat, fields, "read counts")
    if (!nrow(rcmat)) stop("empty read count matrix")
    if (anyNA(rcmat[,fields])) stop("missing input values must be resolved before native fit")
    if (!all(rcmat$Chromosome %in% names(mapping))) stop("unmapped chromosome")
    if (!all(rcmat$catalog_kind %in% c("SNP", "DEPTH")) ||
        !is.logical(rcmat$masked)) stop("explicit catalog kind and mask required")
    for (field in c("Position", "NOR.DP", "NOR.RD", "TUM.DP", "TUM.RD")) {
        x <- rcmat[[field]]
        if (!is.numeric(x) || any(!is.finite(x)) || any(x<0) || any(x!=floor(x)))
            stop("nonnegative integer coordinates/counts required")
    }
    if (any(rcmat$Position<1) || any(rcmat$NOR.RD>rcmat$NOR.DP) ||
        any(rcmat$TUM.RD>rcmat$TUM.DP)) stop("invalid count bounds")
    rcmat$native_chrom <- unname(mapping[rcmat$Chromosome])
    native_rank <- match(rcmat$native_chrom, c(as.character(1:22), "X", "Y"))
    if (is.unsorted(order(native_rank, rcmat$Position)) ||
        anyDuplicated(paste(rcmat$native_chrom, rcmat$Position)))
        stop("read counts must be unique and sorted in native chromosome order")
    rcmat$input_row <- seq_len(nrow(rcmat))
    if (!"catalog_row" %in% names(rcmat)) rcmat$catalog_row <- rcmat$input_row
    if (!is.numeric(rcmat$catalog_row) || any(!is.finite(rcmat$catalog_row)) ||
        any(rcmat$catalog_row<1 | rcmat$catalog_row!=floor(rcmat$catalog_row)) ||
        anyDuplicated(rcmat$catalog_row))
        stop("invalid original catalog row linkage")
    reason <- rep(NA_character_, nrow(rcmat))
    reason[rcmat$native_chrom=="Y"] <- "Y_NOT_SUPPORTED_BY_NATIVE_JOINT_ENGINE"
    reason[rcmat$masked] <- "FROZEN_RESOURCE_MASK"
    reason[is.na(reason) & (rcmat$NOR.DP<config$ndepth |
        rcmat$NOR.DP>=config$ndepthmax)] <- "NORMAL_DEPTH_OUTSIDE_NATIVE_RANGE"
    reason[is.na(reason) & rcmat$TUM.DP==0] <- "NO_SAMPLE_ALLELIC_DEPTH"
    eligible <- rcmat[is.na(reason), , drop=FALSE]
    if (!nrow(eligible)) stop("no eligible rows for joint segmentation")
    eligible$Chromosome <- eligible$native_chrom
    require_fields(bias, c("snp_index", "correction_bias", "bias_unit",
                          "allelic_pseudocount"), "bias")
    if (anyDuplicated(bias$snp_index) || anyNA(bias$snp_index) ||
        anyNA(bias$bias_unit) || !all(bias$bias_unit=="ln_alt_ref") ||
        anyNA(bias$allelic_pseudocount) ||
        !all(abs(bias$allelic_pseudocount-config$allelic_pseudocount)<=
             8*.Machine$double.eps) ||
        !is.numeric(bias$correction_bias) ||
        any(!is.finite(bias$correction_bias) & !is.na(bias$correction_bias)))
        stop("bias units, estimator or locus IDs invalid")
    # pctGCdata 0.3.0 resolves build data through its attached package, not
    # its namespace. Attach only for the real engine, restoring search state.
    if (isNamespace(engine) && identical(unname(getNamespaceName(engine)), "facetsCH") &&
        !"package:pctGCdata" %in% search()) {
        suppressPackageStartupMessages(library("pctGCdata", character.only=TRUE))
        on.exit(detach("package:pctGCdata", unload=FALSE), add=TRUE)
    }
    native <- function(name) get(name, envir=engine, inherits=FALSE)
    set.seed(as.integer(config$seed))
    pmat <- native("procSnps")(eligible, ndepth=config$ndepth,
        het.thresh=config$het_threshold, snp.nbhd=config$snp_neighborhood,
        nX=23, unmatched=TRUE, ndepthmax=config$ndepthmax)
    if (!nrow(pmat)) stop("native SNP filtering produced no rows")
    dmat <- native("counts2logROR")(pmat[pmat$rCountT>0,,drop=FALSE],
                                  gbuild="hg38", unmatched=TRUE)
    key <- paste(dmat$chrom, dmat$maploc, sep="_")
    input_key <- paste(match(eligible$native_chrom,c(as.character(1:22),"X")),
                       eligible$Position,sep="_")
    match_input <- match(key,input_key)
    if (anyNA(match_input) || anyDuplicated(key)) stop("native rows do not map to input loci")
    dmat$valor_uncorrected <- dmat$valor
    dmat$source_input_row <- eligible$input_row[match_input]
    dmat$source_catalog_row <- eligible$catalog_row[match_input]
    dmat$bias_applied <- bias$correction_bias[match(key,bias$snp_index)]
    dmat$valor <- dmat$valor_uncorrected - dmat$bias_applied
    no_allele <- is.na(dmat$valor) | eligible$catalog_kind[match_input]=="DEPTH"
    dmat$valor[no_allele] <- NA_real_
    dmat$het[no_allele] <- 0
    if (any(!is.finite(dmat$cnlr))) stop("nonfinite native logR")
    if (any(!is.finite(dmat$valor) & !is.na(dmat$valor)))
        stop("nonfinite native allelic signal")
    # The native delta scaling divides by mad(diff(logR)); refuse degenerate
    # inputs instead of allowing 0/0. No-allele samples require CNV-only review.
    noise <- mad(diff(dmat$cnlr),na.rm=TRUE)
    if (!is.finite(noise) || noise<=0) stop("insufficient logR dispersion for native scaling")
    if (config$het_scale && sum(dmat$het)==0)
        stop("no informative alleles; requires separate CNV-only review")
    tmp <- native("segsnps")(dmat, cval=config$initial_cval,
        hetscale=config$het_scale, delta=config$delta_cn, gamma=config$gamma)
    prepared <- c(list(pmat=pmat,gbuild="hg38",nX=23),tmp)
    fit <- native("procSample")(prepared, cval=config$final_cval,
                                min.nhet=config$min_nhet)
    retained <- paste(dmat$chrom,dmat$maploc,sep="_")
    eligible_missing <- !input_key %in% retained
    reason[eligible$input_row[eligible_missing]] <- "NATIVE_LOCUS_FILTER"
    locus_audit <- data.frame(input_row=rcmat$input_row,
        catalog_row=rcmat$catalog_row,
        chromosome=rcmat$Chromosome,position_1based=rcmat$Position,
        catalog_kind=rcmat$catalog_kind,normal_depth=rcmat$NOR.DP,
        sample_depth=rcmat$TUM.DP,
        signal_status=ifelse(is.na(reason),"NATIVE_SIGNAL_RETAINED",reason),
        allelic_status="NOT_EVALUATED",stringsAsFactors=FALSE)
    for (i in seq_len(nrow(dmat))) {
        status <- if (eligible$catalog_kind[match_input[i]]=="DEPTH") "DEPTH_ONLY" else
                  if (is.na(dmat$valor_uncorrected[i])) "NATIVE_ALLELIC_SIGNAL_MISSING" else
                  if (is.na(dmat$bias_applied[i])) "BIAS_UNAVAILABLE" else
                  "CORRECTED_ALLELIC_SIGNAL"
        locus_audit$allelic_status[dmat$source_input_row[i]] <- status
    }
    # Export summary is provided separately by facets_segments.R, keeping
    # native segment coordinates and marker spans distinct from evaluability.
    list(status="JOINT_SEGMENTATION_CANDIDATES", fit=fit, signal=dmat,
         locus_audit=locus_audit,
         exclusions=data.frame(input_row=rcmat$input_row[!is.na(reason)],
                               chromosome=rcmat$Chromosome[!is.na(reason)],
                               position_1based=rcmat$Position[!is.na(reason)],
                               reason=reason[!is.na(reason)]),
         config=config, somatic_origin_classified=FALSE,
         copy_neutral_loh_confirmed=FALSE, sensitivity_calibrated=FALSE)
}
