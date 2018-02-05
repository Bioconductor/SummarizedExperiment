### =========================================================================
### collect_rda_objects.R
### -------------------------------------------------------------------------
###
### This script performs STEP 2 of the "Find and update objects" procedure
### described in the README file located in the same folder.
###
### Before you run this script, make sure you performed STEP 1, that is, you
### need to generate input file RDA_FILES. This can be achieved with
### something like:
###
###   cd <dir/you/want/to/search>
###   find . -type d -name '.git' -prune -o -type f -print | \
###       grep -Ei '\.(rda|RData)$' >RDA_FILES
###
### See README file for more information.
###
### Then to run STEP 2 in "batch" mode:
###
###   cd <dir/you/want/to/search>  # RDA_FILES file should be here
###   R CMD BATCH collect_rda_objects.R >collect_rda_objects.log 2>&1 &
###
### This can take a couple of hours to complete...
###
### The output of STEP 2 is a file created in the current directory and named
### RDA_OBJECTS. It has 1 line per serialized object and the 4 following
### fields (separated by tabs):
###   1. Path to rda file (as found in input file RDA_FILES).
###   2. Name of object in rda file.
###   3. Class of object in rda file.
###   4. Package where class of object is defined.
###

INFILE <- "RDA_FILES"
OUTFILE <- "RDA_OBJECTS"

collect_rda_objects <- function(rda_files, outfile="")
{
    cat("", file=outfile)  # create (or overwrite) empty output file
    for (i in seq_along(rda_files)) {
        rda_file <- rda_files[[i]]

        cat("[", i , "/", length(rda_files), "] Loading ", rda_file, " ... ",
            sep="")
        envir <- new.env(parent=emptyenv())
        load(rda_file, envir=envir)
        cat("OK\n")

        for (objname in names(envir)) {
            obj <- get(objname, envir=envir)
            objclass <- class(obj)
            objclass_pkg <- attr(objclass, "package")

            ## Fix 'objclass' and 'objclass_pkg' (they both need to be
            ## character vectors of length 1 before we pass them to paste()).
            objclass <- objclass[[1L]]
            if (is.null(objclass_pkg))
                objclass_pkg <- "."
            outline <- paste(rda_file, objname, objclass, objclass_pkg,
                             sep="\t")
            cat(outline, "\n", sep="", file=OUTFILE, append=TRUE)
        }
    }
}

rda_files <- read.table(INFILE, stringsAsFactors=FALSE)[[1L]]
collect_rda_objects(rda_files, outfile=OUTFILE)

