##
## makeSummarizedExperimentFromExpressionSet

## coercion

.from_rowRanges_to_FeatureData <- function(from)
{
    if (is(from, "GRanges")) {
        fd <- .from_GRanges_to_FeatureData(from)
    } else if (is(from, "GRangesList")) {
        fd <- .from_GRangesList_to_FeatureData(from)
    } else {
        stop("class ", sQuote(class(from)),
             " is not a supported type for rowRanges coercion")
    }
    featureNames(fd) <- names(from)
    fd
}

.from_GRanges_to_FeatureData <- function(from)
{
    data <- as.data.frame(from)

    ## the first mcols are automatically included in the data.frame from
    ## as.data.frame, the secondary mcols holds the metadata for the first
    ## metadata columns.
    metaData <- mcols(mcols(from, use.names=FALSE), use.names=FALSE)
    if (is.null(metaData)) {
        metaData <- as.data.frame(matrix(ncol=0, nrow=NCOL(data)))
    } else {
        metaData <- as.data.frame(metaData)
    }
    AnnotatedDataFrame(data, metaData)
}
.from_GRangesList_to_FeatureData <- function(from)
{
    data <- as.data.frame(mcols(from, use.names=FALSE))

    ## the first mcols are automatically included in the data.frame from
    ## as.data.frame, the secondary mcols holds the metadata for the first
    ## metadata columns.
    metaData <- mcols(mcols(from, use.names=FALSE), use.names=FALSE)
    if (is.null(metaData)) {
        metaData <- as.data.frame(matrix(ncol=0, nrow=NCOL(data)))
    } else {
        metaData <- as.data.frame(metaData)
    }
    AnnotatedDataFrame(data, metaData)
}

.from_AnnotatedDataFrame_to_DataFrame <- function(from)
{
    df <- DataFrame(pData(from), row.names=rownames(from))
    mcols(df) <- DataFrame(varMetadata(from))
    df
}

.from_DataFrame_to_AnnotatedDataFrame <- function(df)
{
    data <- as(df, "data.frame")
    metaData <- mcols(df, use.names=FALSE)
    if (is.null(metaData)) {
        metaData <- as.data.frame(matrix(ncol=0, nrow=NCOL(data)))
    } else {
        metaData <- as(metaData, "data.frame")
    }
    AnnotatedDataFrame(data, metaData)
}

## If the ExpressionSet has featureData with range information make
## GRanges out of that, otherwise make an empty GRangesList with names
## from the featureNames
naiveRangeMapper <- function(from)
{
    nms <- featureNames(from)
    res <- tryCatch({
        makeGRangesFromDataFrame(pData(featureData(from)),
                                 keep.extra.columns = TRUE)
    }, error = function(e) {
        res <- relist(GRanges(), vector("list", length=length(nms)))
        mcols(res) <- .from_AnnotatedDataFrame_to_DataFrame(featureData(from))
        res
    })
    names(res) <- nms
    res
}

# Simple ProbeId to Range mapper
# Probes with multiple ranges are dropped
# The sign of the chromosome location is assumed to contain the strand
# information
probeRangeMapper <- function(from)
{
    annotation <- annotation(from)
    if (identical(annotation, character(0))) {
        return(naiveRangeMapper(from))
    }
    if (requireNamespace("annotate", quietly = TRUE)) {
        annotationPackage <- annotate::annPkgName(annotation)
        test <- require(annotationPackage, character.only = TRUE,
                        quietly = TRUE)
        if (test) {
            db <- get(annotationPackage, envir = asNamespace(annotationPackage))
            pid <- featureNames(from)
            locs <- AnnotationDbi::select(
                db, pid, columns = c("CHR", "CHRLOC", "CHRLOCEND"))
            locs <- na.omit(locs)
            dups <- duplicated(locs$PROBEID)
            if (any(dups)) {
                locs <- locs[!dups, , drop = FALSE]
            }
            strand <- ifelse(locs$CHRLOC > 0, "+", "-")
            res <- GRanges(seqnames = locs$CHR,
                           ranges = IRanges(abs(locs$CHRLOC),
                                            abs(locs$CHRLOCEND)),
                           strand = strand)
            names(res) <- locs$PROBEID

            if (NROW(res) < length(pid)) {
                warning(length(pid) - NROW(res),
                        " probes could not be mapped.", call. = FALSE)
            }
            res
        } else {
            stop("Failed to load ", sQuote(annotationPackage), " package",
                 call. = FALSE)
        }
    } else {
        stop("Failed to load annotate package", call. = FALSE)
    }
}

# Simple ProbeId to Gene mapper
# Is there a way to get the txDb given the annotation package?
geneRangeMapper <- function(txDbPackage, key = "ENTREZID")
{
    function(from) {
        annotation <- annotation(from)
        if (identical(annotation, character(0))) {
            return(naiveRangeMapper(from))
        }
        if (requireNamespace("annotate", quietly = TRUE)) {
            annotationPackage <- annotate::annPkgName(annotation)
            test <- require(annotationPackage, character.only = TRUE,
                            quietly = TRUE)
            if (test) {
                db <- get(annotationPackage,
                          envir = asNamespace(annotationPackage))
                pid <- featureNames(from)
                probeIdToGeneId <-
                    AnnotationDbi::mapIds(db, pid, key, "PROBEID")
                geneIdToProbeId <-
                    setNames(names(probeIdToGeneId), probeIdToGeneId)

                if (requireNamespace(txDbPackage, quietly = TRUE)) {
                    txDb <- get(txDbPackage, envir = asNamespace(txDbPackage))
                    genes <- GenomicFeatures::genes(txDb)
                    probesWithAMatch <-
                        probeIdToGeneId[probeIdToGeneId %in% names(genes)]
                    res <- genes[probesWithAMatch]
                    names(res) <- geneIdToProbeId[names(res)]
                    if (NROW(res) < length(pid)) {
                        warning(length(pid) - NROW(res),
                                " probes could not be mapped.", call. = FALSE)
                    }
                    res
                } else {
                    stop("Failed to load ", sQuote(txDbPackage), " package",
                         call. = FALSE)
                }

            } else {
                stop("Failed to load ", sQuote(annotationPackage), " package",
                     call. = FALSE)
            }
        } else {
            stop("Failed to load annotate package", call. = FALSE)
        }
    }
}

makeSummarizedExperimentFromExpressionSet <-
    function(from, mapFun = naiveRangeMapper, ...)
{
    mapFun <- match.fun(mapFun)
    rowRanges <- mapFun(from, ...)
    matches <- match(names(rowRanges),
                     featureNames(from),
                     nomatch = 0)
    from <- from[matches, drop = FALSE]
    assays <- as.list(assayData(from))
    colData <- .from_AnnotatedDataFrame_to_DataFrame(phenoData(from))
    metadata <- SimpleList(
        experimentData = experimentData(from),
        annotation = annotation(from),
        protocolData = protocolData(from)
    )

    SummarizedExperiment(
        assays = assays,
        rowRanges = rowRanges,
        colData = colData,
        metadata = metadata
    )
}

setAs("ExpressionSet", "RangedSummarizedExperiment", function(from)
{
    makeSummarizedExperimentFromExpressionSet(from)
})

setAs("ExpressionSet", "SummarizedExperiment", function(from)
{
    as(makeSummarizedExperimentFromExpressionSet(from), "SummarizedExperiment")
})

setAs("RangedSummarizedExperiment", "ExpressionSet",
      function(from)
{
    assayData <- list2env(as.list(assays(from)))

    numAssays <- length(assayData)

    if (numAssays == 0) {
        assayData$exprs <- new("matrix")
    } else if (!"exprs" %in% ls(assayData)) {
        ## if there isn't an exprs assay we need to pick one as exprs,
        ## so rename the first element exprs and issue a warning.
        exprs <- ls(assayData)[[1]]
        warning("No assay named ", sQuote("exprs"), " found, renaming ",
                exprs, " to ", sQuote("exprs"), ".")
        assayData[["exprs"]] <- assayData[[exprs]]
        rm(list=exprs, envir=assayData)
    }
    lockEnvironment(assayData, bindings = TRUE)

    featureData <- .from_rowRanges_to_FeatureData(rowRanges(from))
    phenoData <- .from_DataFrame_to_AnnotatedDataFrame(colData(from))

    metadata <- metadata(from)

    experimentData <- if (!is.null(metadata$experimentData)) {
        metadata$experimentData
    } else {
        MIAME()
    }

    annotation <- if (!is.null(metadata$annotation)) {
        metadata$annotation
    } else {
        character()
    }

    protocolData <- if (!is.null(metadata$protocolData)) {
        metadata$protocolData
    } else {
        annotatedDataFrameFrom(assayData, byrow=FALSE)
    }

    ExpressionSet(assayData,
                  phenoData = phenoData,
                  featureData = featureData,
                  experimentData = experimentData,
                  annotation = annotation,
                  protocolData = protocolData
                  )
})

setAs("SummarizedExperiment", "ExpressionSet", function(from)
    as(as(from, "RangedSummarizedExperiment"), "ExpressionSet")
)
