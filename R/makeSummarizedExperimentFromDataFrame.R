### =========================================================================
### makeSummarizedExperimentFromDataFrame()
### -------------------------------------------------------------------------

### 'df' must be a data.frame or DataFrame object.
makeSummarizedExperimentFromDataFrame <-
    function(df,
             ...,
             seqinfo = NULL,
             starts.in.df.are.0based = FALSE)
    {
        rowRanges <- makeGRangesFromDataFrame(
            df,
            ...,
            keep.extra.columns = FALSE,
            seqinfo = seqinfo,
            starts.in.df.are.0based = starts.in.df.are.0based)

        # Find column names for rowRanges
        granges_cols <-
            GenomicRanges:::.find_GRanges_cols(names(df), ...)

        rangedNames <- names(df)[na.omit(granges_cols)]
        idx <- match(rangedNames, names(df))
        counts <- as.matrix(df[, -idx, drop = FALSE])

        if (!is(as.vector(counts), "numeric"))
            stop("failed to coerce non-range columns to 'numeric'")

        SummarizedExperiment(
            assays=SimpleList(counts), rowRanges=rowRanges)
    }
