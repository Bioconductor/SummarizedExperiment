##
## RangedSummarizedExperiment

setClass("RangedSummarizedExperiment",
    contains="Annotated",                     # for metadata slot/getter/setter
    representation(
        rowRanges="GenomicRangesORGRangesList", # rows and their annotations
        colData="DataFrame",                  # columns and their annotations
        assays="Assays"),                     # Data -- e.g., list of matricies
    prototype(
        rowRanges=GRanges(),
        assays=GenomicRanges:::.ShallowSimpleListAssays(data=SimpleList())))

## validity

.valid.SummarizedExperiment.assays_current <- function(x)
{
    if (!is(slot(x, "assays"), "Assays"))
        return("'assays' is out-of-date; use updateObject()")
    NULL
}

.valid.SummarizedExperiment.assays_class <- function(x)
{
    ok <- sapply(assays(x, withDimnames=FALSE), function(cl) {
        (!is.null(dim(cl))) && (length(dim(cl)) >= 2L)
    })
    if (!all(ok))
        return("'assays' must be matrix-like with 2 (or more?) dimensions")
    NULL
}

.valid.RangedSummarizedExperiment.rowRanges_dims <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), nrow) ==
             length(rowRanges(x))))
        return("'rowRanges' length differs from 'assays' nrow")
    NULL
}

.valid.SummarizedExperiment.colData_dims <- function(x)
{
    if (!all(sapply(assays(x, withDimnames=FALSE), ncol) ==
             nrow(colData(x))))
        return("'colData' nrow differs from 'assays' ncol")
    NULL
}

.valid.SummarizedExperiment.assays_dims <- function(x)
{
    c(.valid.RangedSummarizedExperiment.rowRanges_dims(x),
      .valid.SummarizedExperiment.colData_dims(x))
}

.valid.RangedSummarizedExperiment <- function(x)
{
    c(.valid.SummarizedExperiment.assays_current(x),
      msg <- .valid.SummarizedExperiment.assays_class(x),
      if (is.null(msg)) {
          .valid.SummarizedExperiment.assays_dims(x)
      } else NULL)
}

setValidity2("RangedSummarizedExperiment", .valid.RangedSummarizedExperiment)

.GRangesList_assays <-
    function(assays)
{
    m <- assays[[1]]
    n <- nrow(m)
    names <- rownames(m)
    relist(GRanges(), PartitioningByEnd(integer(n), names=names))
}

setMethod(SummarizedExperiment, "SimpleList",
   function(assays, rowRanges=GRangesList(), colData=DataFrame(),
            metadata=list())
{
    if (missing(rowRanges) && 0L != length(assays)) {
        rowRanges <- .GRangesList_assays(assays)
    }
    if (missing(colData) && 0L != length(assays)) {
        nms <- colnames(assays[[1]])
        if (is.null(nms) && 0L != ncol(assays[[1]]))
            stop("'SummarizedExperiment' assay colnames must not be NULL")
        colData <- DataFrame(row.names=nms)
    }

    FUN <- function(x) {
        exp <- list(names(rowRanges), rownames(colData))
        ## dimnames as NULL or list(NULL, NULL)
        all(sapply(dimnames(x), is.null)) ||
            ## or consistent with row / colData
            identical(dimnames(x)[1:2], exp)
    }
    if (!all(sapply(assays, FUN)))
        assays <- endoapply(assays, unname)
    if (!is(assays, "Assays"))
        assays <- GenomicRanges:::.ShallowSimpleListAssays(data=assays)
    new("RangedSummarizedExperiment", metadata=as.list(metadata),
                                      rowRanges=rowRanges,
                                      colData=colData,
                                      assays=assays)
})

setMethod(SummarizedExperiment, "missing",
    function(assays, ...)
{
    SummarizedExperiment(SimpleList(), ...)
})

setMethod(SummarizedExperiment, "list",
    function(assays, ...)
{
    SummarizedExperiment(do.call(SimpleList, assays), ...)
})

setMethod(SummarizedExperiment, "matrix",
    function(assays, ...)
{
    if (is.list(assays))
        ## special case -- matrix of lists
        assays <- list(assays)
    SummarizedExperiment(SimpleList(assays), ...)
})

## update / clone

setMethod(GenomicRanges:::clone, "RangedSummarizedExperiment",  # not exported
    function(x, ...)
{
    ## S4Vectors:::extraArgsAsList would prevent using clone on
    ## subclasses
    args <- list(...)
    firstTime <- TRUE
    for (nm in names(args)) {
        s <- slot(x, nm)
        v <- args[[nm]]
        if (is(s, "ShallowData"))
            v <- GenomicRanges:::clone(s, data=v)
        if (firstTime) {
            slot(x, nm, FALSE) <- v
            firstTime <- FALSE
        } else {
            `slot<-`(x, nm, FALSE, v)
        }
    }
    x
})

setGeneric("value",  # not exported
    function(x, name, ...) standardGeneric("value"),
    signature = "x")

setMethod("value", "RangedSummarizedExperiment",  # not exported
    function(x, name, ...)
{
    s <- slot(x, name)
    if (is(s, "ShallowData"))
        s <- s$data
    s
})

## Simple 'getters' / 'setters'

## We override the default "metadata<-" method (defined in S4Vectors for
## Annotated objects) with a more memory efficient one that uses clone()
## in order to minimize data copy (ShallowData trick).
setReplaceMethod("metadata", "RangedSummarizedExperiment",
    function(x, value)
{
    GenomicRanges:::clone(x, metadata=value)
})

setMethod(rowRanges, "RangedSummarizedExperiment",
    function(x, ...) value(x, "rowRanges"))

.RangedSummarizedExperiment.rowRanges.replace <-
    function(x, ..., value)
{
    x <- GenomicRanges:::clone(x, ..., rowRanges=value)
    msg <- .valid.RangedSummarizedExperiment.rowRanges_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("rowRanges", c("RangedSummarizedExperiment", "GenomicRanges"),
    .RangedSummarizedExperiment.rowRanges.replace)

setReplaceMethod("rowRanges", c("RangedSummarizedExperiment", "GRangesList"),
    .RangedSummarizedExperiment.rowRanges.replace)

setMethod(colData, "RangedSummarizedExperiment",
    function(x, ...) value(x, "colData"))

setReplaceMethod("colData", c("RangedSummarizedExperiment", "DataFrame"),
    function(x, ..., value)
{
    x <- GenomicRanges:::clone(x, ..., colData=value)
    msg <- .valid.SummarizedExperiment.colData_dims(x)
    if (!is.null(msg))
        stop(msg)
    x
})

setMethod("assayNames", "RangedSummarizedExperiment",
    function(x, ...)
{
    names(assays(x, withDimnames=FALSE))
})

setMethod("assayNames<-", c("RangedSummarizedExperiment", "character"),
    function(x, ..., value)
{
    names(assays(x, withDimnames=FALSE)) <- value
    x
})

setMethod(assays, "RangedSummarizedExperiment",
    function(x, ..., withDimnames=TRUE)
{
    if (withDimnames)
        endoapply(value(x, "assays"), "dimnames<-", dimnames(x))
    else
        value(x, "assays")
})

.SummarizedExperiment.assays.replace <-
    function(x, ..., withDimnames=TRUE, value)
{
    ## withDimnames arg allows names(assays(se, withDimnames=FALSE)) <- value
    ok <- vapply(value, function(elt, xdimnames) {
        e <- dimnames(elt)
        (is.null(e[[1]]) || identical(e[[1]], xdimnames[[1]])) &&
             (is.null(e[[2]]) || identical(e[[2]], xdimnames[[2]]))
    }, logical(1), xdimnames=dimnames(x))
    if (!all(ok))
        stop("current and replacement dimnames() differ")
    x <- GenomicRanges:::clone(x, ..., assays=value)
    msg <- .valid.RangedSummarizedExperiment(x)
    if (!is.null(msg))
        stop(msg)
    x
}

setReplaceMethod("assays", c("RangedSummarizedExperiment", "SimpleList"),
    .SummarizedExperiment.assays.replace)

setReplaceMethod("assays", c("RangedSummarizedExperiment", "list"),
    .SummarizedExperiment.assays.replace)

## convenience for common use case
setMethod(assay, c("RangedSummarizedExperiment", "missing"),
    function(x, i, ...)
{
    assays <- assays(x, ...)
    if (0L == length(assays))
        stop("'assay(<", class(x), ">, i=\"missing\", ...) ",
             "length(assays(<", class(x), ">)) is 0'")
    assays[[1]]
})

setMethod(assay, c("RangedSummarizedExperiment", "numeric"),
    function(x, i, ...)
{
    tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop("'assay(<", class(x), ">, i=\"numeric\", ...)' ",
             "invalid subscript 'i'\n", conditionMessage(err))
    })
})

setMethod(assay, c("RangedSummarizedExperiment", "character"),
    function(x, i = names(x)[1], ...)
{
    msg <- paste0("'assay(<", class(x), ">, i=\"character\", ...)' ",
                  "invalid subscript 'i'")
    res <- tryCatch({
        assays(x, ...)[[i]]
    }, error=function(err) {
        stop(msg, "\n", conditionMessage(err))
    })
    if (is.null(res))
        stop(msg, "\n'i' not in names(assays(<", class(x), ">))")
    res
})

setReplaceMethod("assay", c("RangedSummarizedExperiment", "missing", "matrix"),
    function(x, i, ..., value)
{
    if (0L == length(assays(x)))
        stop("'assay(<", class(x), ">) <- value' ", "length(assays(<",
             class(x), ">)) is 0")
    assays(x)[[1]] <- value
    x
})

setReplaceMethod("assay",
    c("RangedSummarizedExperiment", "numeric", "matrix"),
    function(x, i = 1, ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

setReplaceMethod("assay",
    c("RangedSummarizedExperiment", "character", "matrix"),
    function(x, i = names(x)[1], ..., value)
{
    assays(x, ...)[[i]] <- value
    x
})

## cannonical location for dim, dimnames; dimnames should be checked
## for consistency (if non-null) and stripped from assays on
## construction, or added from assays if row/col names are NULL in
## <RangedSummarizedExperiment> but not assays. dimnames need to be added on
## to assays when assays() invoked
setMethod(dim, "RangedSummarizedExperiment",
    function(x)
{
    c(length(rowRanges(x)), nrow(colData(x)))
})

setMethod(dimnames, "RangedSummarizedExperiment",
    function(x)
{
    list(names(rowRanges(x)), rownames(colData(x)))
})

setReplaceMethod("dimnames", c("RangedSummarizedExperiment", "list"),
    function(x, value)
{
    rowRanges <- rowRanges(x)
    names(rowRanges) <- value[[1]]
    colData <- colData(x)
    rownames(colData) <- value[[2]]
    GenomicRanges:::clone(x, rowRanges=rowRanges, colData=colData)
})

setReplaceMethod("dimnames", c("RangedSummarizedExperiment", "NULL"),
    function(x, value)
{
    dimnames(x) <- list(NULL, NULL)
    x
})

## Subset -- array-like

.SummarizedExperiment.charbound <-
    function(idx, txt, fmt)
{
    orig <- idx
    idx <- match(idx, txt)
    if (any(bad <- is.na(idx))) {
        msg <- paste(BiocGenerics:::selectSome(orig[bad]), collapse=" ")
        stop(sprintf(fmt, msg))
    }
    idx
}

.SummarizedExperiment.assays.subset <-
    function(x, i, j, ...)
{
    ## need to expand Rle's for subsetting standard matrix
    if (!missing(i) && !missing(j)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[i, j, drop=FALSE],
                   x[i, j, , drop=FALSE],
                   x[i, j, , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    } else if (!missing(i)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[i, , drop=FALSE],
                   x[i, , , drop=FALSE],
                   x[i, , , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    } else if (!missing(j)) {
        fun <- function(x) {
            switch(length(dim(x)),
                   stop("'[' on assays() with 1 dimension not supported"),
                   x[, j, drop=FALSE],
                   x[, j, , drop=FALSE],
                   x[, j, , , drop=FALSE],
                   stop("'[' on assays() with >4 dimensions not supported"))
        }
    }
    endoapply(assays(x, withDimnames=FALSE), fun)
}

setMethod("[", c("RangedSummarizedExperiment", "ANY", "ANY"),
    function(x, i, j, ..., drop=TRUE)
{
    if (1L != length(drop) || (!missing(drop) && drop))
        warning("'drop' ignored '[,", class(x), ",ANY,ANY-method'")

    if (missing(i) && missing(j))
        return(x)

    if (!missing(i) && is.character(i)) {
        fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
        i <- .SummarizedExperiment.charbound(i, rownames(x), fmt)
    }
    if (!missing(j) && is.character(j)) {
        fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
        j <- .SummarizedExperiment.charbound(j, colnames(x), fmt)
    }

    if (!missing(i) && !missing(j)) {
        ii <- as.vector(i)
        jj <- as.vector(j)
        x <- GenomicRanges:::clone(x, ..., rowRanges=rowRanges(x)[i],
            colData=colData(x)[j, , drop=FALSE],
            assays=.SummarizedExperiment.assays.subset(x, ii, jj))
    } else if (missing(i)) {
        jj <- as.vector(j)
        x <- GenomicRanges:::clone(x, ..., colData=colData(x)[j, , drop=FALSE],
            assays=.SummarizedExperiment.assays.subset(x, j=jj))
    } else {                            # missing(j)
        ii <- as.vector(i)
        x <- GenomicRanges:::clone(x, ..., rowRanges=rowRanges(x)[i],
            assays=.SummarizedExperiment.assays.subset(x, ii))
    }
    x
})

.SummarizedExperiment.assays.subsetgets <-
    function(x, i, j, ..., value)
{
    ## need to expand Rle's for subsetting standard matrix
    if (!missing(i) && !missing(j)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[i, j] <- value,
                   x[i, j, ] <- value,
                   x[i, j, , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    } else if (!missing(i)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[i, ] <- value,
                   x[i, , ] <- value,
                   x[i, , , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    } else if (!missing(j)) {
        fun <- function(x, value) {
            switch(length(dim(x)),
                   stop("'[<-' on assays() with 1 dimension not supported"),
                   x[, j] <- value,
                   x[, j, ] <- value,
                   x[, j, , ] <- value,
                   stop("'[<-' on assays() with >4 dimensions not supported"))
            x
        }
    }
    a <- assays(x, withDimnames=FALSE)
    v <- assays(value, withDimnames=FALSE)
    mendoapply(fun, x=a, value=v)
}

setReplaceMethod("[",
    c("RangedSummarizedExperiment", "ANY", "ANY", "RangedSummarizedExperiment"),
    function(x, i, j, ..., value)
{
    if (missing(i) && missing(j))
        return(value)

    if (!missing(i) && is.character(i)) {
        fmt <- paste0("<", class(x), ">[i,] index out of bounds: %s")
        i <- .SummarizedExperiment.charbound(i, rownames(x), fmt)
    }

    if (!missing(j) && is.character(j)) {
        fmt <- paste0("<", class(x), ">[,j] index out of bounds: %s")
        j <- .SummarizedExperiment.charbound(j, colnames(x), fmt)
    }

    if (!missing(i) && !missing(j)) {
        ii <- as.vector(i)
        jj <- as.vector(j)
        x <- GenomicRanges:::clone(x, ..., metadata=c(metadata(x),
                                                      metadata(value)),
            rowRanges=local({
                r <- rowRanges(x)
                r[i] <- rowRanges(value)
                names(r)[ii] <- names(rowRanges(value))
                r
            }), colData=local({
                c <- colData(x)
                c[j,] <- colData(value)
                rownames(c)[jj] <- rownames(colData(value))
                c
            }), assays=.SummarizedExperiment.assays.subsetgets(x, ii, jj,
                  ..., value=value))
        msg <- .valid.SummarizedExperiment.assays_dims(x)
    } else if (missing(i)) {
        jj <- as.vector(j)
        x <- GenomicRanges:::clone(x, ..., metadata=c(metadata(x),
                                                      metadata(value)),
            colData=local({
                c <- colData(x)
                c[j,] <- colData(value)
                rownames(c)[jj] <- rownames(colData(value))
                c
            }), assays=.SummarizedExperiment.assays.subsetgets(x, j=jj,
                  ..., value=value))
        msg <- .valid.SummarizedExperiment.colData_dims(x)
    } else {                            # missing(j)
        ii <- as.vector(i)
        x <- GenomicRanges:::clone(x, ..., metadata=c(metadata(x),
                                                      metadata(value)),
            rowRanges=local({
                r <- rowRanges(x)
                r[i] <- rowRanges(value)
                names(r)[ii] <- names(rowRanges(value))
                r
            }), assays=.SummarizedExperiment.assays.subsetgets(x, ii,
                  ..., value=value))
        msg <- .valid.RangedSummarizedExperiment.rowRanges_dims(x)
    }
    if (!is.null(msg))
        stop(msg)
    x
})

## rbind, cbind

## Appropriate for objects with different ranges and same samples.
setMethod("rbind", "RangedSummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .rbind.RangedSummarizedExperiment(args)
})

.rbind.RangedSummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, colnames)))
            stop("'...' objects must have the same colnames")
    if (!.compare(lapply(args, ncol)))
            stop("'...' objects must have the same number of samples")
    rowRanges <- do.call(c, lapply(args, rowRanges))
    colData <- .cbind.DataFrame(args, colData, "colData")
    assays <- .bind.arrays(args, rbind, "assays")
    metadata <- do.call(c, lapply(args, metadata))

    initialize(args[[1]], assays=assays, rowRanges=rowRanges,
               colData=colData, metadata=metadata)
}

## Appropriate for objects with same ranges and different samples.
setMethod("cbind", "RangedSummarizedExperiment",
    function(..., deparse.level=1)
{
    args <- unname(list(...))
    .cbind.RangedSummarizedExperiment(args)
})

.cbind.RangedSummarizedExperiment <- function(args)
{
    if (!.compare(lapply(args, rowRanges), TRUE))
        stop("'...' object ranges (rows) are not compatible")
    rowRanges <- rowRanges(args[[1]])
    mcols(rowRanges) <- .cbind.DataFrame(args, mcols, "mcols")
    colData <- do.call(rbind, lapply(args, colData))
    assays <- .bind.arrays(args, cbind, "assays")
    metadata <- do.call(c, lapply(args, metadata))

    initialize(args[[1]], assays=assays, rowRanges=rowRanges,
               colData=colData, metadata=metadata)
}

.compare <- function(x, GenomicRanges=FALSE)
{
    x1 <- x[[1]]
    if (GenomicRanges) {
        if (is(x1, "GRangesList")) {
            x <- lapply(x, unlist)
            x1 <- x[[1]]
        }
        for (i in seq_along(x)[-1]) {
            if (length(x1) != length(x[[i]]))
                return(FALSE)
            ok <- x1 == x[[i]]
            if (!all(ok))
                return(FALSE)
        }
        return(TRUE)
    } else {
        all(sapply(x[-1],
            function(xelt) all(identical(xelt, x[[1]]))))
    }
}

.cbind.DataFrame <- function(args, accessor, accessorName)
{
    lst <- lapply(args, accessor)
    if (!.compare(lst)) {
        nms <- lapply(lst, names)
        nmsv <- unlist(nms, use.names=FALSE)
        names(nmsv) <- rep(seq_along(nms), elementLengths(nms))
        dups <- duplicated(nmsv)
        ## no duplicates
        if (!any(dups))
            return(do.call(cbind, lst))
        ## confirm duplicates are the same
        lapply(nmsv[duplicated(nmsv)], function(d) {
            if (!.compare(lapply(lst, "[", d)))
                stop("column(s) '", unname(d),
                     "' in ", sQuote(accessorName),
                     " are duplicated and the data do not match")})
        ## remove duplicates
        do.call(cbind, lst)[,!dups]
    } else {
        lst[[1]]
    }
}

.bind.array.elements <- function(index, lst, bind) {
    e1 <- lapply(lst, "[[", index)
    dim <- .get.assay.dimension(e1, bind)
    if (is.na(dim[3])) {
        do.call(bind, e1)
    } else {
        e2 <- lapply(1:dim[3], function(i) {
            do.call(bind, lapply(e1, "[", ,,i))
        })
        array(do.call(c, e2), dim=dim)
    }
}

.bind.arrays <- function(args, bind, accessor)
{
    lst <- lapply(args, accessor)
    if (!length(elts <- unique(elementLengths(lst))) == 1L)
        stop("elements in ", sQuote(accessor),
             " must have the same length")
    if (elts == 0L)
        return(GenomicRanges:::.ShallowSimpleListAssays(data=SimpleList()))
    var <- lapply(lst,  names)
    if (is.null(uvar <- unique(unlist(var)))) {
        ## no names, match by position
        res <- lapply(seq_along(elts), .bind.array.elements, lst=lst, bind=bind)
    } else {
        ## match by name
        if (!.compare(var))
            stop("elements in ", sQuote(accessor),
                 " must have the same names")
        res <- lapply(uvar, .bind.array.elements, lst=lst, bind=bind)
        names(res) <- uvar
    }
    GenomicRanges:::.ShallowSimpleListAssays(data=SimpleList(res))
}

.get.assay.dimension <- function(lst, bind)
{
    dim <- lapply(lst, dim)
    if (!.compare(lapply(dim, "[", 3)))
        stop("elements in assays must have the same dimension")
    if (identical(bind, cbind))
        c(dim[[1]][1], do.call(sum, lapply(dim, "[", 2)), dim[[1]][3])
    else
        c(do.call(sum, lapply(dim, "[", 1)), dim[[1]][2], dim[[1]][3])
}

## $, $<-, [[, [[<- for colData access

setMethod("[[", c("RangedSummarizedExperiment", "ANY", "missing"),
    function(x, i, j, ...)
{
    colData(x)[[i, ...]]
})

setReplaceMethod("[[",
    c("RangedSummarizedExperiment", "ANY", "missing", "ANY"),
    function(x, i, j, ..., value)
{
    colData(x)[[i, ...]] <- value
    x
})

.DollarNames.RangedSummarizedExperiment <- function(x, pattern)
    grep(pattern, names(colData(x)), value=TRUE)

setMethod("$", "RangedSummarizedExperiment",
    function(x, name)
{
    colData(x)[[name]]
})

setReplaceMethod("$", c("RangedSummarizedExperiment", "ANY"),
    function(x, name, value)
{
    colData(x)[[name]] <- value
    x
})

## show

setMethod(show, "RangedSummarizedExperiment",
    function(object)
{
    selectSome <- BiocGenerics:::selectSome
    scat <- function(fmt, vals=character(), exdent=2, ...)
    {
        vals <- ifelse(nzchar(vals), vals, "''")
        lbls <- paste(BiocGenerics:::selectSome(vals), collapse=" ")
        txt <- sprintf(fmt, length(vals), lbls)
        cat(strwrap(txt, exdent=exdent, ...), sep="\n")
    }
    cat("class:", class(object), "\n")
    cat("dim:", dim(object), "\n")
    expt <- names(metadata(object))
    if (is.null(expt))
        expt <- character(length(metadata(object)))
    scat("metadata(%d): %s\n", expt)
    nms <- names(assays(object, withDimnames=FALSE))
    if (is.null(nms))
        nms <- character(length(assays(object, withDimnames=FALSE)))
    scat("assays(%d): %s\n", nms)
    dimnames <- dimnames(object)
    dlen <- sapply(dimnames, length)
    if (dlen[[1]]) scat("rownames(%d): %s\n", dimnames[[1]])
    else scat("rownames: NULL\n")
    scat("rowRanges metadata column names(%d): %s\n",
         names(mcols(rowRanges(object))))
    if (dlen[[2]]) scat("colnames(%d): %s\n", dimnames[[2]])
    else cat("colnames: NULL\n")
    scat("colData names(%d): %s\n", names(colData(object)))
})

## compatibility

.from_SummarizedExperiment_to_RangedSummarizedExperiment <- function(from)
    new("RangedSummarizedExperiment", metadata=as.list(from@exptData),
                                      rowRanges=from@rowData,
                                      colData=from@colData,
                                      assays=from@assays)

### Should work on any object that (1) belongs to a class that now derives
### from RangedSummarizedExperiment and (2) was created at a time when it was
### deriving from the old SummarizedExperiment class defined in GenomicRanges.
### For example: DESeqDataSet objects created before the definition of class
### DESeqDataSet was modified to extend RangedSummarizedExperiment instead of
### SummarizedExperiment.  
setMethod(updateObject, "RangedSummarizedExperiment",
    function(object, ..., verbose=FALSE)
{
    has_SummarizedExperiment_internal_structure <-
      all(sapply(slotNames("SummarizedExperiment"), .hasSlot, object=object))
    if (!has_SummarizedExperiment_internal_structure)
        return(object)
    rse <- .from_SummarizedExperiment_to_RangedSummarizedExperiment(object)
    xslotnames <- setdiff(slotNames(class(object)), slotNames(class(rse)))
    xslots <- attributes(object)[xslotnames]
    #do.call("new", c(list(Class=class(object), rse), xslots))

    ## The line above doesn't work because of a bug in R (see
    ## https://stat.ethz.ch/pipermail/r-devel/2015-May/071130.html),
    ## so we use the workaround below.
    rse_slots <- attributes(rse)[slotNames(class(rse))]
    do.call("new", c(list(Class=class(object)), rse_slots, xslots))
})

