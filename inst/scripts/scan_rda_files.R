### =========================================================================
### scan_rda_files.R
### -------------------------------------------------------------------------
###
### This script performs STEP 2 of the "Find and update objects" procedure
### described in the 'Find_and_update_objects.txt' file located in the same
### folder.
###
### Before you run this script, make sure you performed STEP 1, that is, you
### need to generate input file 'rda_files'. This can be achieved with
### something like:
###
###   cd <dir/you/want/to/search>
###   find . -type d -name '.svn' -prune -o -type f -print | \
###       grep -Ei '\.(rda|RData)$' >rda_files
###
### See Find_and_update_objects.txt for more information.
###
### Then to run STEP 2 in "batch" mode:
###
###   cd <dir/you/want/to/search>  # the 'rda_files' file should be here
###   R CMD BATCH scan_rda_files.R >scan_rda_files.log 2>&1 &
###
### This can take a couple of hours to complete...
###
### The output of STEP 2 is a file created in the current directory and named
### 'rda_objects'. It has 1 line per serialized object and the 4 following
### fields (separated by tabs):
###   1. Path to rda file (as found in input file 'rda_files').
###   2. Name of object in rda file.
###   3. Class of object in rda file.
###   4. Package where class of object is defined.
###

INFILE <- "rda_files"
OUTFILE <- "rda_objects"

scanRdaFiles <- function(rda_files, outfile="")
{
    cat("", file=outfile)  # create (or overwrite) empty output file
    for (i in seq_along(rda_files)) {
        rda_path <- rda_files[[i]]

        cat("[", i , "/", length(rda_files), "] Loading ", rda_path, " ... ",
            sep="")
        envir <- new.env(parent=emptyenv())
        load(rda_path, envir=envir)
        cat("OK\n")

        for (objname in names(envir)) {
            obj <- get(objname, envir=envir)
            obj_class <- class(obj)
            obj_class_pkg <- attr(obj_class, "package")

            ## Fix 'obj_class' and 'obj_class_pkg' (they both need to be
            ## character vectors of length 1 before we pass them to paste()).
            obj_class <- obj_class[[1L]]
            if (is.null(obj_class_pkg))
                obj_class_pkg <- "."
            outline <- paste(rda_path, objname, obj_class, obj_class_pkg,
                             sep="\t")
            cat(outline, "\n", sep="", file=OUTFILE, append=TRUE)
        }
    }
}

rda_files <- read.table(INFILE, stringsAsFactors=FALSE)[[1L]]
scanRdaFiles(rda_files, outfile=OUTFILE)

