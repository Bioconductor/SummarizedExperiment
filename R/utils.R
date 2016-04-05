### =========================================================================
### Some low-level utils
### -------------------------------------------------------------------------
###
### Unless stated otherwise, nothing in this file is exported.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Bind arrays with an arbitrary number of dimensions along the 1st or 2nd
### dimension
###

### Return a matrix with one row per dim and one column per object.
.get_and_check_objects_dims <- function(objects, no.check.along)
{
    dims <- lapply(objects, dim)
    ndims <- lengths(dims)
    ndim <- ndims[[1L]]
    if (!all(ndims == ndim))
        stop("all the objects to bind must have the same nb of dimensions")
    tmp <- unlist(dims, use.names=FALSE)
    if (is.null(tmp))
        stop("the objects to bind have no dimensions")
    dims <- matrix(tmp, nrow=ndim)
    tmp <- dims[-no.check.along, , drop=FALSE]
    if (!all(tmp == tmp[ , 1L]))
        stop("objects to bind have incompatible dimensions")
    dims
}

### 'objects' is assumed to be a list of vector-like objects.
### 'block_lens' is assumed to be an integer vector parallel to 'objects'
### specifying the block length for each object in 'objects'. In addition the
### length of 'object[[i]]' must be 'k * block_lens[[i]]' (k is the same for
### all the objects).
.intertwine_blocks <- function(objects, block_lens)
{
    nobject <- length(objects)
    objects_lens <- lengths(objects)
    objects_cumlens <- cumsum(objects_lens)
    k <- objects_lens[[1L]] %/% block_lens[[1L]]
    x <- unlist(objects, recursive=FALSE, use.names=FALSE)
    ranges <- lapply(seq_len(nobject),
        function(i) {
            width <- block_lens[[i]]
            offset <- if (i == 1L) 0L else objects_cumlens[[i - 1L]]
            successiveIRanges(rep.int(width, k), from=offset + 1L)
        })
    ranges <- do.call(c, ranges)
    i <- as.vector(matrix(seq_len(nobject * k), nrow=nobject, byrow=TRUE))
    extractROWS(x, ranges[i])
}

.combine_dimnames <- function(objects, dims, along)
{
    ndim <- nrow(dims)
    dimnames <- lapply(seq_len(ndim),
        function(n) {
            for (x in objects) {
                dn <- dimnames(x)[[n]]
                if (!is.null(dn))
                    return(dn)
            }
            NULL
        })
    along_names <- lapply(objects, function(x) dimnames(x)[[along]])
    along_names_lens <- lengths(along_names)
    if (any(along_names_lens != 0L)) {
        fix_idx <- which(along_names_lens != dims[along, ])
        along_names[fix_idx] <- lapply(dims[along, fix_idx], character)
    }
    along_names <- unlist(along_names, use.names=FALSE)
    if (!is.null(along_names))
        dimnames[[along]] <- along_names
    if (all(S4Vectors:::sapply_isNULL(dimnames)))
        dimnames <- NULL
    dimnames
}

### A stripped-down version of abind::abind().
### Some differences:
###   (a) Treatment of dimnames: .simple_abind() treatment of dimnames is
###       consistent with base::rbind() and base::cbind(). This is not the
###       case for abind::abind() which does some strange things with the
###       dimnames.
###   (b) Performance: .simple_abind() has a little bit more overhead than
###       abind::abind(). This makes it slower on small objects. However it
###       tends to be slightly faster on big objects.
.simple_abind <- function(..., along)
{
    objects <- list(...)
    object_is_NULL <- S4Vectors:::sapply_isNULL(objects)
    if (any(object_is_NULL))
        objects <- objects[!object_is_NULL]
    if (length(objects) == 0L)
        return(NULL)
    if (length(objects) == 1L)
        return(objects[[1L]])

    ## Check dim compatibility.
    dims <- .get_and_check_objects_dims(objects, no.check.along=along)

    ## Perform the binding.
    block_lens <- dims[along, ]
    for (n in seq_len(along - 1L))
        block_lens <- block_lens * dims[n, ]
    ans <- .intertwine_blocks(objects, block_lens)

    ## Set the dim.
    ans_dim <- dims[ , 1L]
    ans_dim[[along]] <- sum(dims[along, ])
    dim(ans) <- ans_dim

    ## Combine and set the dimnames.
    dimnames(ans) <- .combine_dimnames(objects, dims, along=along)
    ans
}

arbind_default <- function(...) .simple_abind(..., along=1L)
acbind_default <- function(...) .simple_abind(..., along=2L)

