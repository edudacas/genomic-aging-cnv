# File boundary only. The Python orchestrator must verify stage provenance first.
# TSV configuration is parsed as data, never source()/dget()/eval().
read_exact_table <- function(path, columns, classes) {
    x <- read.delim(path, header=TRUE, sep="\t", quote="", comment.char="",
        check.names=FALSE, stringsAsFactors=FALSE, colClasses=classes,
        na.strings="NA", fill=FALSE)
    if (!identical(names(x),columns)) stop("unexpected TSV schema")
    x
}

read_joint_files <- function(matrix_path,bias_path,config_path,map_path) {
    counts <- read_exact_table(matrix_path,
        c("Chromosome","Position","NOR.DP","NOR.RD","TUM.DP","TUM.RD",
          "catalog_kind","masked","catalog_row"),
        c("character",rep("numeric",5),"character","logical","numeric"))
    bias <- read_exact_table(bias_path,
        c("snp_index","correction_bias","bias_unit","allelic_pseudocount"),
        c("character","numeric","character","numeric"))
    raw <- read_exact_table(config_path,c("key","value"),c("character","character"))
    numeric_fields <- c("ndepth","ndepthmax","het_threshold","snp_neighborhood",
        "initial_cval","final_cval","gamma","min_nhet","seed","allelic_pseudocount",
        "sample_het_min_depth","delta_cn")
    if (anyNA(raw) || anyDuplicated(raw$key) ||
        !setequal(raw$key,c(numeric_fields,"reference_build","het_scale")))
        stop("explicit exact joint configuration required")
    config <- as.list(setNames(raw$value,raw$key))
    for (field in numeric_fields) {
        value <- suppressWarnings(as.numeric(config[[field]]))
        if (length(value)!=1 || !is.finite(value)) stop("non-numeric config value")
        config[[field]] <- value
    }
    if (!config$het_scale %in% c("TRUE","FALSE")) stop("explicit boolean required")
    config$het_scale <- config$het_scale=="TRUE"
    mapping <- read_exact_table(map_path,c("source","native"),c("character","character"))
    if (!nrow(mapping) || anyNA(mapping) || anyDuplicated(mapping$source) ||
        anyDuplicated(mapping$native)) stop("invalid chromosome mapping")
    config$chromosome_map <- setNames(mapping$native,mapping$source)
    list(counts=counts,bias=bias,config=config)
}

run_joint_files <- function(matrix_path,bias_path,config_path,map_path,output,
                            engine=asNamespace("facetsCH")) {
    inputs <- read_joint_files(matrix_path,bias_path,config_path,map_path)
    if (file.exists(output)) stop("existing attempt must be preserved")
    if (!dir.create(output,recursive=TRUE)) stop("cannot create joint output")
    writeLines("JOINT_FIT_STARTED",file.path(output,"STARTED.txt"))
    tryCatch({
        result <- run_grch38_joint(inputs$counts,inputs$bias,inputs$config,engine)
        segments <- summarize_native_segments(result$fit$jointseg,inputs$config$chromosome_map)
        write_tsv <- function(x,name) write.table(x,file.path(output,name),sep="\t",
            row.names=FALSE,col.names=TRUE,quote=FALSE,na="NA")
        write_tsv(segments,"segments.tsv")
        write_tsv(result$locus_audit,"locus_audit.tsv")
        write_tsv(result$exclusions,"exclusions.tsv")
        write_tsv(result$signal,"signal.tsv")
        saveRDS(result,file.path(output,"native_fit.rds"))
        writeLines(capture.output(sessionInfo()),file.path(output,"R_session.txt"))
        # Completion is last. Outer orchestration hashes these files and must
        # never treat this native text marker as a clinical/validated result.
        writeLines(c("JOINT_SEGMENTATION_CANDIDATES","sensitivity_calibrated=FALSE",
            "somatic_origin_classified=FALSE","copy_neutral_loh_confirmed=FALSE"),
            file.path(output,"COMPLETE.txt"))
        invisible(result)
    },error=function(e) {
        writeLines(conditionMessage(e),file.path(output,"FAILURE.txt"))
        stop(e)
    })
}
