# Contains methods for combineRows and combineCols. These serve as more
# fault-tolerant relaxed counterparts to rbind and cbind, respectively.

setMethod("combineRows", "SummarizedExperiment", function(x, y, ..., use.names=TRUE, delayed=TRUE) {
    all.se <- list(x, y, ...)

    # Combining the rowData.
    all.rd <- lapply(all.se, rowData) 
    tryCatch({
        com.rd <- do.call(combineRows, all.rd)
    }, error=function(e) {
        stop(paste0("failed to combine rowData of SummarizedExperiment objects:\n  ", conditionMessage(e)))
    })

    # Combining the colData. This constructs mappings of the columns for each
    # SE to the columns of the final object.
    all.cd <- lapply(all.se, colData)
    tryCatch({
        com.cd <- do.call(combineCols, c(all.cd, list(use.names=use.names)))
    }, error=function(e) {
        stop(paste0("failed to combine colData of SummarizedExperiment objects:\n  ", conditionMessage(e)))
    })

    if (use.names) {
        # If combineCols succeeded, all SE's should have valid row names.
        all.names <- rownames(com.cd)
        mappings <- vector("list", length(all.se))
        for (i in seq_along(mappings)) {
            mappings[[i]] <- match(all.names, rownames(all.cd[[i]]))
        }
    } else {
        mappings <- NULL
    }

    # Assembling the SE.
    SummarizedExperiment(
        assays=.unsplit_assays(all.se, mappings, delayed=delayed),
        colData=com.cd,
        rowData=com.rd,
        metadata=unlist(lapply(all.se, metadata), recursive=FALSE, use.names=FALSE)
    )
})

.unsplit_assays <- function(all.se, mappings, delayed) {
    all.assays <- lapply(all.se, assays, withDimnames=FALSE)
    each.assay.names <- lapply(all.assays, names)
    no.assay.names <- vapply(each.assay.names, is.null, TRUE)
    
    create_dummy <- function(nr, nc) {
        if (!delayed) {
            array(c(nr, nc), data=NA)
        } else {
            ConstantArray(c(nr, nc), value=NA)
        }
    }

    inflate_existing <- function(mat, i) {
        if (delayed) {
            mat <- DelayedArray(mat)
        }
        if (!is.null(i)) {
            absent <- is.na(i)
            if (any(absent)) {
                i[absent] <- ncol(mat)+1L
                mat <- cbind(mat, create_dummy(nrow(mat), 1L))
            }
            mat <- mat[,i,drop=FALSE]
        }
        mat
    }

    if (any(no.assay.names)) {
        if (!all(no.assay.names)) {
            stop("named and unnamed assays cannot be mixed")
        }
        n.assays <- unique(lengths(all.assays))
        if (length(n.assays)!=1L) {
            stop("all SummarizedExperiments should have the same number of unnamed assays")
        }
        for (s in seq_along(all.se)) {
            cur.assays <- all.assays[[s]]
            i <- mappings[[s]]

            # Filling in all missing columns.
            for (a in seq_len(n.assays)) {
                mat <- inflate_existing(cur.assays[[a]], i)
                cur.assays[[a]] <- mat
            }
            all.assays[[s]] <- cur.assays
        }
    } else {
        all.assay.names <- Reduce(union, each.assay.names)
        for (s in seq_along(all.se)) {
            cur.se <- all.se[[s]]
            cur.assays <- all.assays[[s]]
            i <- mappings[[s]]

            # Filling in all missing assay names and columns.
            for (a in all.assay.names) {
                if (a %in% names(cur.assays)) {
                    mat <- inflate_existing(cur.assays[[a]], i)
                } else {
                    if (is.null(i)) {
                        nc <- ncol(cur.se)
                    } else {
                        nc <- length(i)
                    }
                    mat <- create_dummy(nrow(cur.se), nc)
                }
                cur.assays[[a]] <- mat
            }
            all.assays[[s]] <- cur.assays
        }
    }

    # Re-use assay rbind'ing machinery.
    all.assays <- lapply(all.assays, Assays)
    combined <- do.call(rbind, all.assays)
    as(combined, "SimpleList") 
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

    # Filling in the empties. We promote everyone to a GRL if any of the
    # individual entries are GRLs.
    has.grl <- any(vapply(final.rd, FUN=is, class2="GRangesList", FUN.VALUE=TRUE))
    for (e in which(vapply(final.rd, is.null, FUN.VALUE=TRUE))) {
        cur.se <- all.se[[e]]
        if (has.grl) {
            empty <- GRangesList(rep(list(GRanges()), nrow(cur.se)))
        } else {
            # TODO: a better convention for missing intervals?
            empty <- GRanges(rep("unknown:1-0", nrow(cur.se))) 
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
