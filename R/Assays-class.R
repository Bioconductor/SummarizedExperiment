### =========================================================================
### Assays objects
### -------------------------------------------------------------------------
###
### The Assays API consists of:
###   (a) The Assays() constructor function.
###   (b) Lossless back and forth coercion from/to SimpleList. The coercion
###       method from SimpleList doesn't need (and should not) validate the
###       returned object.
###   (c) The following methods, split in 2 groups:
###       - List-like methods:   length, names, names<-, getListElement, and
###                              setListElement
###       - Matrix-like methods: dim, [, [<-, rbind, cbind
###
### An Assays concrete subclass needs to implement (b) (required) plus
### optionally any of the methods in (c). The reason they are optionals
### is that default methods are provided and they work on any Assays
### derivative as long as lossless back and forth coercion from/to
### SimpleList works.
###
### IMPORTANT
### ---------
###
### 1. Nobody in the Assays hierarchy is allowed to inherit from SimpleList
###    because of the conflicting semantic of [.
###
### 2. Methods that return a modified Assays object (a.k.a. endomorphisms),
###    that is, [ as well as replacement methods names<-, setListElement,
###    and [<-, must respect the copy-on-change contract. With objects that
###    don't make use of references internally, the developer doesn't need
###    to take any special action for that because it's automatically taken
###    care of by R itself. However, for objects that do make use of
###    references internally (e.g. environments, external pointers, pointer
###    to a file on disk, etc...), the developer needs to be careful to
###    implement endomorphisms with copy-on-change semantics. This can
###    easily be achieved (and is what the default methods for Assays objects
###    do) by performaing a full (deep) copy of the object before modifying it
###    instead of trying to modify it in-place. However note that this full
###    (deep) copy can be very expensive and is actually not necessary in
###    order to achieve copy-on-change semantics: it's enough (and often
###    preferrable for performance reasons) to copy only the parts of the
###    object that need to be modified.


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Assays class
###

setClass("Assays", contains="RectangularData", representation("VIRTUAL"))

### Validity

.valid.Assays <- function(x)
{
    assays <- try(as(x, "SimpleList"), silent=TRUE)
    if (inherits(assays, "try-error"))
        return("'as(x, \"SimpleList\")' must work")
    if (!is(assays, "SimpleList"))
        return("'as(x, \"SimpleList\")' must return a SimpleList object")
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

### Always return a SimpleList object by default. Will return a NULL only
### if 'as.null.if.no.assay' is set to TRUE and no assays are supplied.
normarg_assays <- function(assays, as.null.if.no.assay=FALSE)
{
    if (!isTRUEorFALSE(as.null.if.no.assay))
        stop(wmsg("'as.null.if.no.assay' must be TRUE or FALSE"))

    if (is.null(assays)) {
        if (as.null.if.no.assay)
            return(NULL)
        return(SimpleList())
    }

    ## The truth is that the assays can be any array-like objects
    ## with at least 2 dimensions, not just matrix-like objects.
    error_msg <- c("'assays' must be a list or SimpleList of ",
                   "matrix-like elements, or a matrix-like object, ",
                   "or a NULL (see '?SummarizedExperiment')")

    if (is(assays, "Assays"))
        stop(wmsg(error_msg))

    assays_dim <- dim(assays)
    ## Some objects like SplitDataFrameList have a dim() method that
    ## returns a non-NULL object (a matrix!) even though they don't have
    ## an array-like semantic.
    if (!is.matrix(assays_dim) && length(assays_dim) >= 2L)
        #return(SimpleList(assays))  # broken on a data frame
        return(new2("SimpleList", listData=list(assays), check=FALSE))

    if (!is(assays, "SimpleList")) {
        if (is.list(assays)) {
            #assays <- do.call(SimpleList, assays) # broken on a list of
                                                   # data frames
            assays <- new2("SimpleList", listData=assays, check=FALSE)
        } else if (is(assays, "List")) {
            assays <- as(assays, "SimpleList")  # could fail
        } else {
            stop(wmsg(error_msg))
        }
    }
    if (length(assays) == 0L && as.null.if.no.assay)
        return(NULL)
    assays
}

### Always return a SimpleAssays object by default. Will return a NULL only
### if 'as.null.if.no.assay' is set to TRUE and no assays are supplied.
Assays <- function(assays=SimpleList(), as.null.if.no.assay=FALSE)
{
    if (!isTRUEorFALSE(as.null.if.no.assay))
        stop(wmsg("'as.null.if.no.assay' must be TRUE or FALSE"))
    ## Starting with SummarizedExperiment 1.15.4, we wrap the user-supplied
    ## assays in a SimpleAssays object instead of a ShallowSimpleListAssays
    ## object. Note that there are probably hundreds (if not thousands) of
    ## serialized SummarizedExperiment objects around that use
    ## ShallowSimpleListAssays. These objects should keep working as before!
    if (!is(assays, "SimpleAssays")) {
        if (is(assays, "Assays")) {
            ## Will turn any Assays derivative (e.g. ShallowSimpleListAssays)
            ## into a SimpleAssays object.
            assays <- as(as(assays, "SimpleList"), "SimpleAssays")
        } else {
            assays <- normarg_assays(assays, as.null.if.no.assay)
            if (is.null(assays))
                return(NULL)
            assays <- as(assays, "SimpleAssays")
            validObject(assays)
        }
    }
    if (length(assays) == 0L && as.null.if.no.assay)
        return(NULL)
    assays
}

### updateObject

.updateObject_Assays <- function(object, ..., verbose=FALSE)
{
    assays <- as(object, "SimpleList")
    assays <- endoapply(assays,
        function(assay)
            updateObject(assay, ..., verbose=verbose)
    )
    if (length(assays) == 0L)
        return(NULL)
    as(assays, "SimpleAssays")
}

setMethod("updateObject", "Assays", .updateObject_Assays)

### Accessors

setMethod("length", "Assays",
    function(x)
    {
        x <- as(x, "SimpleList")
        callGeneric()
    }
)

setMethod("names", "Assays",
    function(x)
    {
        x <- as(x, "SimpleList")
        callGeneric()
    }
)

setReplaceMethod("names", "Assays",
    function(x, value)
    {
        ans_class <- class(x)
        x <- as(x, "SimpleList")
        as(callGeneric(), ans_class)
    }
)

setMethod("getListElement", "Assays",
    function(x, i, exact=TRUE)
    {
        x <- as(x, "SimpleList")
        callGeneric()
    }
)

setMethod("setListElement", "Assays",
    function(x, i, value)
    {
        ans_class <- class(x)
        x <- as(x, "SimpleList")
        ans <- as(callGeneric(), ans_class)
        validObject(ans)
        ans
    }
)

setMethod("dim", "Assays",
    function(x)
    {
        if (length(x) == 0L)
            return(c(0L, 0L))
        dim(getListElement(x, 1L))[1:2]
    }
)

### 2D-Subsetting

### Subset each assay in Assays object 'x' along its first 2 dimensions.
### If not missing, 'i' and/or 'j' are assumed to be valid subscripts that
### can be used to subset an ordinary matrix or array.
.extract_Assays_subset <- function(x, i, j)
{
    subscripts12 <- list(if (missing(i)) quote(expr=) else i,
                         if (missing(j)) quote(expr=) else j)
    extract_assay_subset <- function(a) {
        ndim <- length(dim(a))
        stopifnot(ndim >= 2L)  # should never happen
        more_subscripts <- rep.int(list(quote(expr=)), ndim - 2L)
        args <- c(list(a), subscripts12, more_subscripts, list(drop=FALSE))
        do.call(`[`, args)
    }
    assays <- as(x, "SimpleList")
    as(endoapply(assays, extract_assay_subset), class(x))
}

setMethod("[", "Assays",
    function(x, i, j, ..., drop=TRUE) .extract_Assays_subset(x, i, j)
)

### Subassign each assay in Assays object 'x' along its first 2 dimensions.
### If not missing, 'i' and/or 'j' are assumed to be valid subscripts that
### can be used to subassign an ordinary matrix or array.
.replace_Assays_subset <- function(x, i, j, value)
{
    subscripts12 <- list(if (missing(i)) quote(expr=) else i,
                         if (missing(j)) quote(expr=) else j)
    replace_assay_subset <- function(a, v) {
        ndim <- length(dim(a))
        stopifnot(ndim >= 2L)  # should never happen
        more_subscripts <- rep.int(list(quote(expr=)), ndim - 2L)
        args <- c(list(a), subscripts12, more_subscripts, list(value=v))
        do.call(`[<-`, args)
    }
    assays <- as(x, "SimpleList")
    values <- as(value, "SimpleList")
    as(mendoapply(replace_assay_subset, assays, values), class(x))
}

setReplaceMethod("[", "Assays",
    function(x, i, j, ..., value) .replace_Assays_subset(x, i, j, value)
)

### rbind/cbind

### 'assays' is assumed to be an unnamed list of length >= 1
.bind_assays <- function(assays, along.cols=FALSE)
{
    if (length(dim(getListElement(assays, 1L))) == 2L) {
        BINDING_FUN <- if (along.cols) "cbind" else "rbind"
    } else {
        BINDING_FUN <- if (along.cols) "acbind" else "arbind"
    }
    do.call(BINDING_FUN, assays)
}

.bind_Assays_objects <- function(objects, along.cols=FALSE)
{
    if (length(objects) == 0L)
        return(Assays())
    lens <- sapply(objects, length)
    if (length(unique(lens)) != 1)
        stop("the objects to bind must have the same number of assays")
    len1 <- lens[1L]
    if (len1 == 0L)
        return(Assays())
    var <- lapply(objects, names)
    uvar <- unique(unlist(var))
    if (is.null(uvar)) {
        ## no names, match by position
        res <- lapply(seq_len(len1), function(index) {
            assays <- lapply(objects, getListElement, index)
            .bind_assays(assays, along.cols=along.cols)
        })
    } else {
        ## match by name
        ok <- all(vapply(var, function(x, y) identical(sort(x), y),
                         logical(1), sort(uvar)))
        if (!ok)
            stop("assays must have the same names()")
        res <- lapply(uvar, function(index) {
            assays <- lapply(objects, getListElement, index)
            .bind_assays(assays, along.cols=along.cols)
        })
        names(res) <- uvar
    }
    as(SimpleList(res), class(getListElement(objects, 1L)))
}

setMethod("rbind", "Assays",
    function(..., deparse.level=1)
    {
        objects <- unname(list(...))
        .bind_Assays_objects(objects, along.cols=FALSE)
    }
)

setMethod("cbind", "Assays",
    function(..., deparse.level=1)
    {
        objects <- unname(list(...))
        .bind_Assays_objects(objects, along.cols=TRUE)
    }
)

### Having "arbind" and "acbind" methods for Matrix objects will make rbind()
### and cbind() work on Assays objects with Matrix list elements.
### Maybe these methods should be defined next to the arbind() and acbind()
### generics (which are defined in the IRanges package) but that would require
### to make IRanges depend on the Matrix package.
setMethod("arbind", "Matrix", function(...) rbind(...))
setMethod("acbind", "Matrix", function(...) cbind(...))


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### SimpleAssays class
###
### Store the Assays in a SimpleList object.
###

### SimpleAssays cannot contain SimpleList because of the conflicting
### semantic of [.
setClass("SimpleAssays",
    contains="Assays",
    representation(data="SimpleList")
)

### We only need to implement the REQUIRED coercions.

setAs("SimpleList", "SimpleAssays",
    function(from) new2("SimpleAssays", data=from, check=FALSE)
)

setAs("SimpleAssays", "SimpleList", function(from) from@data)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### ShallowSimpleListAssays class
###
### WARNING: Looks like reference classes as implemented in the methods
### package are a bit problematic e.g. all.equal() can return false negatives
### after a serialization/deserialization cycle on a ref object as reported
### here https://stat.ethz.ch/pipermail/bioc-devel/2019-May/015112.html
### Anyway their use in the assays slot of a SummarizedExperiment object is
### probably not needed anymore now that R is "shallow copy by default"
### according to this comment by Michael:
###     https://github.com/Bioconductor/SummarizedExperiment/issues/16#issuecomment-455541415
###

.ShallowData <- setRefClass("ShallowData",
    fields = list( data = "ANY" ))

.ShallowSimpleListAssays0 <- setRefClass("ShallowSimpleListAssays",
    fields = list( data = "SimpleList" ),
    contains = c("ShallowData", "Assays"))

### We only need to implement the REQUIRED coercions.

setAs("SimpleList", "ShallowSimpleListAssays",
    function(from) .ShallowSimpleListAssays0(data=from)
)

setAs("ShallowSimpleListAssays", "SimpleList", function(from) from$data)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### AssaysInEnv class
###
### A *broken* alternative to ShallowSimpleListAssays that does NOT respect
### the copy-on-change contract (only provided for illustration purposes).
###
### We implement the REQUIRED coercions plus OPTIONAL methods: length, names,
### names<-, getListElement, and setListElement.
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
        value <- S4Vectors:::normarg_names(value, class(x), length(x))
        x@envir[[.NAMES_SYMBOL]] <- value
        x
    }
)

setMethod("getListElement", "AssaysInEnv",
    function(x, i, exact=TRUE)
    {
        key <- setNames(ls(x@envir, sorted=TRUE), names(x))[[i]]
        get(key, envir=x@envir)
    }
)

### Does NOT respect the copy-on-change contract!
setMethod("setListElement", "AssaysInEnv",
    function(x, i, value)
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

