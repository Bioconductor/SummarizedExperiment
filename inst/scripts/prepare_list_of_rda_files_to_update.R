### =========================================================================
### prepare_list_of_rda_files_to_update.R
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
###   cd <dir/you/want/to/scan>
###   find . -type d -name '.svn' -prune -o -type f -print | \
###       grep -Ei '\.(rda|RData)$' >rda_files
###
### See Find_and_update_objects.txt for more information.
###
### Then to run STEP 2 in "batch" mode:
###
###   cd <dir/you/want/to/scan>  # the 'rda_files' file should be here
###   R CMD BATCH prepare_list_of_rda_files_to_update.R \
###              >prepare_list_of_rda_files_to_update.log 2>&1 &
###
### This can take a few hours to complete...
###
### The output of STEP 2 is a file created in the current directory and named
### 'rda_files_to_update'. It has 1 line per object to update and 4 fields
### separated by tabs:
###   1. Path to rda file (as found in input file 'rda_files').
###   2. Name of object to update.
###   3. Class of object to update.
###   4. Package where class of object to update is defined.
###

INFILE <- "rda_files"
OUTFILE <- "rda_files_to_update"

library(BiocInstaller)
available_pkgs <- rownames(available.packages(contrib.url(biocinstallRepos())))

installAndLoadPkg <- function(package)
{
    if (suppressWarnings(suppressPackageStartupMessages(
            require(package, character.only=TRUE, quietly=TRUE)
        )))
        return()
    suppressWarnings(suppressPackageStartupMessages(
        library(BiocInstaller, quietly=TRUE)
    ))
    biocLite(package)
    suppressWarnings(suppressPackageStartupMessages(
        library(package, character.only=TRUE, quietly=TRUE)
    ))
}

library(SummarizedExperiment)

prepareListOfRdaFilesToUpdate <- function(rda_files, outfile="")
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
            class_pkg <- attr(class(obj), "package")
            if (is.null(class_pkg) || !(class_pkg %in% available_pkgs))
                next
            installAndLoadPkg(class_pkg)
            if (is(obj, "SummarizedExperiment")
             || is(obj, "RangedSummarizedExperiment")) {
                cat(rda_path, "\t", objname, "\t", class(obj), "\t",
                    class_pkg, "\n", sep="", file=outfile, append=TRUE)
            }
        }
    }
}

rda_files <- read.table(INFILE, stringsAsFactors=FALSE)[[1L]]
prepareListOfRdaFilesToUpdate(rda_files, outfile=OUTFILE)

