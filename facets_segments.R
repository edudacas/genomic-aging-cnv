# Explicit native-to-BED segment summaries. These are evidence, not calls.
summarize_native_segments <- function(signal, chromosome_map) {
    if (!is.character(chromosome_map) || is.null(names(chromosome_map)) ||
        anyNA(chromosome_map) || anyNA(names(chromosome_map)) ||
        any(names(chromosome_map)=="") || anyDuplicated(names(chromosome_map)) ||
        anyDuplicated(unname(chromosome_map)))
        stop("one-to-one explicit chromosome map required")
    required <- c("chrom","maploc","seg","cnlr","valor","het")
    if (!all(required %in% names(signal)) || !nrow(signal))
        stop("nonempty native joint signal required")
    for (field in c("chrom","maploc","seg","cnlr","het")) {
        if (!is.numeric(signal[[field]]) || any(!is.finite(signal[[field]])))
            stop("invalid native segment field")
    }
    for (field in c("chrom","maploc","seg"))
        if (any(signal[[field]]<1 | signal[[field]]!=floor(signal[[field]])))
            stop("integer native coordinate/segment required")
    if (!all(signal$het %in% c(0,1)) ||
        any(!is.finite(signal$valor) & !is.na(signal$valor)))
        stop("invalid allelic evidence")
    if (anyDuplicated(paste(signal$chrom,signal$maploc)) ||
        is.unsorted(order(signal$chrom,signal$maploc)))
        stop("duplicate/unsorted native marker")
    native_names <- c(as.character(1:22),"X","Y")
    source_chrom <- names(chromosome_map)[match(native_names[signal$chrom],chromosome_map)]
    if (length(source_chrom)!=nrow(signal) || anyNA(source_chrom))
        stop("native chromosome lacks explicit reverse mapping")
    group <- paste(signal$chrom,signal$seg,sep="_")
    if (anyDuplicated(rle(group)$values)) stop("discontinuous segment ID")
    rows <- lapply(unique(group),function(key) {
        idx <- which(group==key)
        x <- signal[idx,,drop=FALSE]
        allele <- x$valor[x$het==1 & !is.na(x$valor)]
        start <- min(x$maploc)-1
        end <- max(x$maploc)
        data.frame(chromosome=source_chrom[idx[1]],native_segment=x$seg[1],
            native_start_1based=min(x$maploc),native_end_1based=max(x$maploc),
            start_0based=start,end_0based_exclusive=end,
            marker_bounding_span_bp=end-start,n_markers=nrow(x),
            n_informative_allelic=length(allele),
            mean_log2_ratio=mean(x$cnlr),median_log2_ratio=median(x$cnlr),
            mean_abs_ln_allelic=if(length(allele)) mean(abs(allele)) else NA_real_,
            evidence_class=if(length(allele)) "JOINT_SIGNAL" else "CNV_ONLY_SIGNAL",
            event_class="UNCLASSIFIED",evaluable_bp=NA_real_,
            somatic_origin_classified=FALSE,stringsAsFactors=FALSE)
    })
    do.call(rbind,rows)
}
