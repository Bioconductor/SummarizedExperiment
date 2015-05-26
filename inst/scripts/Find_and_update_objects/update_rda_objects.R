### =========================================================================
### update_rda_objects.R
### -------------------------------------------------------------------------
###
### This script performs STEP 4 of the "Find and update objects" procedure
### described in the README file located in the same folder.
###
### Before you run this script, make sure you performed STEPS 1 & 2 & 3.
### See README file for more information.
###
### Then to run STEP 4 in "batch" mode:
###
###   cd <dir/you/want/to/search>  # RDA_OBJECTS_TO_UPDATE file should be here
###   R CMD BATCH update_rda_objects.R >update_rda_objects.log 2>&1 &
###

INFILE <- "RDA_OBJECTS_TO_UPDATE"

library(SummarizedExperiment)

.update_objects <- function(envir)
{
    updated <- FALSE
    for (objname in names(envir)) {
        obj <- get(objname, envir=envir, inherits=FALSE)
        objclass <- class(obj)
        objclass_pkg <- attr(objclass, "package")
        if (!is.null(objclass_pkg)) {
            suppressWarnings(suppressPackageStartupMessages(
                library(objclass_pkg, character.only=TRUE, quietly=TRUE)
            ))
        }
        if (!SummarizedExperiment:::.has_SummarizedExperiment_internal_structure(obj))
            next()

        cat("  Updating ", objname, " ... ", sep="")
        obj <- updateObject(obj)
        validObject(obj, complete=TRUE)
        assign(objname, obj, envir=envir, inherits=FALSE)
        cat("OK\n")
        updated <- TRUE
    }
    updated
}

updateRdaObjects <- function(rda_objects)
{
    rda_files <- unique(rda_objects[ , "rda_file"])
    for (i in seq_along(rda_files)) {
        rda_file <- rda_files[[i]]

        cat("[", i , "/", length(rda_files), "] Loading ", rda_file, " ... ",
            sep="")
        envir <- new.env(parent=emptyenv())
        load(rda_file, envir=envir)
        cat("OK\n")

        if (.update_objects(envir)) {
            cat("  Saving updated objects to ", rda_file, " ... ", sep="")
            save(list=names(envir), file=rda_file, envir=envir, compress="xz")
            cat("OK\n")
        }
    }
}

rda_objects_to_update <- read.table(INFILE, stringsAsFactors=FALSE)
colnames(rda_objects_to_update) <- c("rda_file", "objname",
                                     "objclass", "objclass_pkg")
updateRdaObjects(rda_objects_to_update)

