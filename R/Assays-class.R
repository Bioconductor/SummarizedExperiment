### =========================================================================
### Assays objects
### -------------------------------------------------------------------------
###
### The Assays, ShallowData, and ShallowSimpleListAssays classes are helper
### classes that are currently defined in the SummarizedExperiment package.
### Move them here once the "old" SummarizedExperiment class in GenomicRanges
### is defunct (in BioC 2.3).
###
### The Assays API consists of:
###   (a) The Assays() constructor function.
###   (b) Lossless back and forth coercion from/to SimpleList. The coercion
###       method from SimpleList doesn't need (and should not) validate the
###       returned object.
###   (c) length, names, names<-, [[, [[<-, dim, [, [<-, rbind, cbind
###
### An Assays croncrete subclass needs to implement (b) (required) plus
### optionally any of the methods in (c).
###
### IMPORTANT: Methods that return a modified Assays object (a.k.a.
### endomorphisms), that is, [ as well as replacement methods names<-, [[<-,
### and [<-, must respect the copy-on-change contract. With objects that
### don't make use of references internally, the developer doesn't need to
### take any special action for that because it's automatically taken care of
### by R itself. However, for objects that do make use of references internally
### (e.g. environments, external pointers, pointer to a file on disk, etc...),
### the developer needs to be careful to implement endomorphisms with
### copy-on-change semantics. This can easily be achieved (and is what the
### default methods for Assays objects do) by performaing a full (deep) copy
### of the object before modifying it instead of trying to modify it in-place.
### Note that the full (deep) copy is not always necessary in order to achieve
### copy-on-change semantics: it's enough (and often preferrable for
### performance reasons) to copy only the parts of the objects that need to
### be modified.
###


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Assays class
###

### Validity

.valid.Assays <- function(x)
{
    assays <- as(x, "SimpleList", strict=FALSE)
    if (!is(assays, "SimpleList"))
        return("'assays' must be a SimpleList object")
    if (length(assays) == 0L)
        return(NULL)

    ## Check dims.
    all_dims <- sapply(assays, function(assay) dim(assay)[1:2])
    if (any(is.na(all_dims)))
        return(wmsg("all assays must be matrix-like objects ",
                    "with 2 (or more?) dimensions"))
    if (!all(all_dims == all_dims[ , 1L]))
        stop("all assays must have the same nrow and ncol")

    NULL
}

setValidity2("Assays", .valid.Assays)

### Constructor

.normarg_assays <- function(assays)
{
    if (!is(assays, "SimpleList")) {
        if (is.list(assays) || is.array(assays)) {
            assays <- SimpleList(assays)
        } else {
            stop("'assays' must be a SimpleList, list or array")
        }
    }
    assays
}

Assays <- function(assays=SimpleList())
{
    assays <- .normarg_assays(assays)
    #ans <- as(assays, "SimpleListAssays")
    ans <- as(assays, "ShallowSimpleListAssays")
    #ans <- as(assays, "AssaysInEnv")  # a *broken* alternative
    validObject(ans)
    ans
}

### Accessors

.SL_get_length <- selectMethod("length", "SimpleList")
setMethod("length", "Assays",
    function(x)
    {
        assays <- as(x, "SimpleList", strict=FALSE)
        .SL_get_length(assays)
    }
)

.SL_get_names <- selectMethod("names", "SimpleList")
setMethod("names", "Assays",
    function(x)
    {
        assays <- as(x, "SimpleList", strict=FALSE)
        .SL_get_names(assays)
    }
)

.SL_set_names <- selectMethod("names<-", "SimpleList")
setReplaceMethod("names", "Assays",
    function(x, value)
    {
        assays <- as(x, "SimpleList", strict=FALSE)
        assays <- .SL_set_names(assays, value)
        as(assays, class(x))
    }
)

setMethod("[[", "Assays",
    function(x, i, j, ...)
    {
        assays <- as(x, "SimpleList", strict=FALSE)
        getListElement(assays, i)
    }
)

setReplaceMethod("[[", "Assays",
    function(x, i, j, ..., value)
    {
        assays <- as(x, "SimpleList", strict=FALSE)
        assays <- setListElement(assays, i, value)
        ans <- as(assays, class(x))
        validObject(ans)
        ans
    }
)

setMethod("dim", "Assays",
    function(x)
    {
        if (length(x) == 0L)
            return(c(0L, 0L))
        dim(x[[1L]])
    }
)

### 2D-Subsetting

.extract_Assays_subset <- function(x, i, j)
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
    assays <- as(x, "SimpleList", strict=FALSE)
    as(endoapply(assays, fun), class(x))
}

setMethod("[", "Assays",
    function(x, i, j, ..., drop=TRUE) .extract_Assays_subset(x, i, j)
)

.replace_Assays_subset <- function(x, i, j, value)
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
    a <- as(x, "SimpleList", strict=FALSE)
    v <- as(value, "SimpleList", strict=FALSE)
    as(mendoapply(fun, x=a, value=v), class(x))
}

setReplaceMethod("[", "Assays",
    function(x, i, j, ..., value) .replace_Assays_subset(x, i, j, value)
)

### rbind/cbind

.get_assay_dimension <- function(lst, bind)
{
    dim <- lapply(lst, dim)
    if (length(unique(sapply(dim, "[", 3))) != 1L)
        stop("elements in assays must have the same dimension")
    if (identical(bind, cbind))
        c(dim[[1]][1], do.call(sum, lapply(dim, "[", 2)), dim[[1]][3])
    else
        c(do.call(sum, lapply(dim, "[", 1)), dim[[1]][2], dim[[1]][3])
}

.bind_assay_elements <- function(index, lst, bind) {
    e1 <- lapply(lst, "[[", index)
    dim <- .get_assay_dimension(e1, bind)
    if (is.na(dim[3])) {
        do.call(bind, e1)
    } else {
        e2 <- lapply(1:dim[3], function(i) {
            do.call(bind, lapply(e1, "[", ,,i))
        })
        array(do.call(c, e2), dim=dim)
    }
}

.bind_Assays <- function(lst, bind)
{
    if (length(lst) == 0L)
        return(Assays())
    lens <- sapply(lst, length)
    len1 <- lens[1L]
    if (any(lens != len1))
        stop("elements in assays must have the same length")
    if (len1 == 0L)
        return(Assays())
    var <- lapply(lst, names)
    if (is.null(uvar <- unique(unlist(var)))) {
        ## no names, match by position
        res <- lapply(seq_along(len1),
                      .bind_assay_elements, lst=lst, bind=bind)
    } else {
        ## match by name
        ok <- all(sapply(var[-1],
                  function(xelt) all(identical(xelt, var[[1]]))))
        if (!ok)
            stop("elements in assays must have the same names")
        res <- lapply(uvar, .bind_assay_elements, lst=lst, bind=bind)
        names(res) <- uvar
    }
    as(res, class(len1))
}

setMethod("rbind", "Assays",
    function(..., deparse.level=1) .bind_Assays(unname(list(...)), rbind)
)

setMethod("cbind", "Assays",
    function(..., deparse.level=1) .bind_Assays(unname(list(...)), cbind)
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SimpleListAssays class
###

### The order of inheritance is important: first Assays, then SimpleList!
setClass("SimpleListAssays", contains=c("Assays", "SimpleList"))

### Lossless back and forth coercion from/to SimpleList are automatically
### taken care of by automatic methods defined by the methods package.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ShallowSimpleListAssays class
###
### We implement the REQUIRED coercions only.
###

### GenomicRanges:::.ShallowSimpleListAssays0 seems a little bit faster than
### two-step constructor GenomicRanges:::.ShallowSimpleListAssays (June 2015)
setAs("SimpleList", "ShallowSimpleListAssays",
    function(from) GenomicRanges:::.ShallowSimpleListAssays0(data=from)
)

setAs("ShallowSimpleListAssays", "SimpleList", function(from) from$data)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### AssaysInEnv class
###
### A *broken* alternative to ShallowSimpleListAssays that does NOT respect
### the copy-on-change contract (only provided for illustration purposes).
###
### We implement the REQUIRED coercions plus OPTIONAL methods: length, names,
### names<-, [[, and [[<-.
###

setClass("AssaysInEnv",
    contains="Assays",
    representation(envir="environment")
)

.NAMES_SYMBOL <- ".names"  # must begin with a . so is ommitted by ls() 

setMethod("length", "AssaysInEnv", function(x) length(x@envir) - 1L)

setMethod("names", "AssaysInEnv", function(x) x@envir[[.NAMES_SYMBOL]])

### Does NOT respect the copy-on-change contract!
setReplaceMethod("names", "AssaysInEnv",
    function(x, value)
    {
        value <- S4Vectors:::normalize_names_replacement_value(value, x)
        x@envir[[.NAMES_SYMBOL]] <- value
        x
    }
)

setMethod("[[", "AssaysInEnv",
    function(x, i, j, ...)
    {
        key <- setNames(ls(x@envir, sorted=TRUE), names(x))[[i]]
        get(key, envir=x@envir)
    }
)

### Does NOT respect the copy-on-change contract!
setReplaceMethod("[[", "AssaysInEnv",
    function(x, i, j, ..., value)
    {
        key <- setNames(ls(x@envir, sorted=TRUE), names(x))[[i]]
        assign(key, value, envir=x@envir)
        x
    }
)

setAs("SimpleList", "AssaysInEnv",
    function(from)
    {
        from <- as.list(from)
        from_names <- names(from)
        keys <- paste(sprintf("%09d", seq_along(from)), from_names, sep=":")
        names(from) <- keys
        envir <- list2env(from, parent=emptyenv())
        envir[[.NAMES_SYMBOL]] <- from_names
        new("AssaysInEnv", envir=envir)
    }
)

setAs("AssaysInEnv", "SimpleList",
    function(from)
        SimpleList(setNames(as.list(from@envir, sorted=TRUE), names(from)))
)

