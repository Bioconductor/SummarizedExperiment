# Contains methods for combineRows and combineCols. These serve as more
# fault-tolerant relaxed counterparts to rbind and cbind, respectively.

setMethod("combineRows", "SummarizedExperiment", function(x, y, ..., use.names=TRUE, delayed=TRUE) {
    all.se <- list(x, y, ...)

    # Combining the rowData.
    all.rd <- lapply(all.se, rowData) 
    all.rd <- do.call(rbind, all.rd) # TODO: replace with combineRows for a DF.

    # Combining the colData. This constructs mappings of the columns for each
    # SE to the columns of the final object.
    all.cd <- lapply(all.se, colData)

    if (use.names) {
        col.names <- lapply(all.se, colnames)
        all.names <- Reduce(union, col.names)

        mappings <- vector("list", length(all.se))
        for (i in seq_along(mappings)) {
            cur.names <- col.names[[i]]
            if (is.null(cur.names) || anyDuplicated(cur.names)) {
                stop(wmsg("each SummarizedExperiment must have non-NULL, ",
                          "non-duplicated column names when 'use.names=TRUE'"))
            }
            mappings[[i]] <- match(all.names, cur.names)
        }

        for (i in seq_along(all.cd)) {
            expanded <- all.cd[[i]][mappings[[i]],,drop=FALSE]
            rownames(expanded) <- all.names
            all.cd[[i]] <- expanded
        }
        all.cd <- do.call(cbind, all.cd) 

    } else {
        nc <- unique(vapply(all.se, ncol, 0L))
        if (length(nc)!=1) {
            stop(wmsg("all SummarizedExperiment objects must have the same ",
                      "number of columns when 'use.names=FALSE'"))
        }
        mappings <- NULL

        all.cd <- do.call(cbind, all.cd) 
    }

    # Assembling the SE.
    SummarizedExperiment(
        assays=.unsplit_assays(all.se, mappings, delayed=delayed),
        colData=all.cd,
        rowData=all.rd,
        metadata=unlist(lapply(all.se, metadata), recursive=FALSE, use.names=FALSE)
    )
})

.unsplit_assays <- function(all.se, mappings, delayed) {
    all.assays <- unique(unlist(lapply(all.se, assayNames)))
    combined <- vector("list", length(all.assays))
    names(combined) <- all.assays

    for (a in all.assays) {
        current <- vector("list", length(all.se))

        for (s in seq_along(all.se)) {
            cur.se <- all.se[[s]]
            i <- mappings[[s]]

            if (a %in% assayNames(cur.se)) {
                mat <- assay(cur.se, a, withDimnames=FALSE)
                if (delayed) {
                    mat <- DelayedArray(mat)
                }

                # Doing something to handle the missing columns, if there are any.
                if (!is.null(i)) {
                    absent <- is.na(i)
                    if (any(absent)) {
                        i[absent] <- ncol(mat)+1L
                        mat <- cbind(mat, ConstantArray(c(nrow(cur.se), 1L), value=NA))
                    }
                    mat <- mat[,i,drop=FALSE]
                }
            } else {
                if (is.null(i)) {
                    nc <- ncol(cur.se)
                } else {
                    nc <- length(i)
                }

                mat <- ConstantArray(c(nrow(cur.se), nc), value=NA)
                if (!delayed) {
                    mat <- as.matrix(mat)
                }
            }

            current[[s]] <- mat
        }
        combined[[a]] <- do.call(rbind, current)
    }

    combined
}

setMethod("combineRows", "RangedSummarizedExperiment", function(x, y, ..., use.names=TRUE, delayed=TRUE) {
    out <- callNextMethod()

    all.se <- list(x, y, ...)
    final.rd <- vector("list", length(all.se))
    for (s in seq_along(all.se)) {
        cur.se <- all.se[[s]]
        rr <- rowRanges(cur.se)
        mcols(rr) <- NULL
        final.rd[[s]] <- rr
    }

    # Filling in the empties. We promote everyone to a GRL if any of the individuals are GRLs.
    has.grl <- any(vapply(final.rd, FUN=is, class2="GRangesList", FUN.VALUE=TRUE))
    for (e in which(vapply(final.rd, is.null, FUN.VALUE=TRUE))) {
        cur.se <- all.se[[e]]
        if (has.grl) {
            empty <- GRangesList(rep(list(GRanges()), nrow(cur.se)))
        } else {
            empty <- GRanges(rep("unknown:1-0", nrow(cur.se))) # TODO: a better convention for missing intervals?
        }
        mcols(empty) <- rowData(cur.se)
        names(empty) <- rownames(cur.se)
        final.rd[[e]] <- empty
    }

    rr <- do.call(c, final.rd)
    mcols(rr) <- rowData(out)
    rowRanges(out) <- rr
    out 
})
