source("facets_segments.R")
x <- data.frame(chrom=c(1,1,1,23),maploc=c(1,100,200,5),seg=c(1,2,2,3),
                cnlr=c(-0.4,0,0.02,0),valor=c(0.2,NA,NA,0),het=c(1,0,0,1))
y <- summarize_native_segments(x,c(chr1="1",chrX="X"))
stopifnot(nrow(y)==3,y$start_0based[1]==0,y$end_0based_exclusive[1]==1,
    y$marker_bounding_span_bp[1]==1,y$marker_bounding_span_bp[2]==101,
    y$n_markers[2]==2,is.na(y$mean_abs_ln_allelic[2]),
    y$evidence_class[2]=="CNV_ONLY_SIGNAL",all(is.na(y$evaluable_bp)),
    all(y$event_class=="UNCLASSIFIED"),y$chromosome[3]=="chrX")
bad <- x; bad$maploc[3] <- 100
stopifnot(inherits(try(summarize_native_segments(bad,c(chr1="1",chrX="X")),silent=TRUE),"try-error"))
bad <- x; bad$seg[1:3] <- c(1,2,1)
stopifnot(inherits(try(summarize_native_segments(bad,c(chr1="1",chrX="X")),silent=TRUE),"try-error"))
cat("NATIVE_SEGMENT_EXPORT_CONTRACT_PASS\n")
