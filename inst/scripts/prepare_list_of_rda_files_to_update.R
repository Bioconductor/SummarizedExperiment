### =========================================================================
### prepare_list_of_rda_files_to_update.R
### -------------------------------------------------------------------------
###
### This script performs STEP 3 of the "Find and update objects" procedure
### described in the 'Find_and_update_objects.txt' file located in the same
### folder.
###
### Before you run this script, make sure you performed STEP 1 & 2.
### See Find_and_update_objects.txt for more information.
###
### Then to run STEP 3 in "batch" mode:
###
###   cd <dir/you/want/to/search>  # the 'rda_files' file should be here
###   R CMD BATCH prepare_list_of_rda_files_to_update.R \
###              >prepare_list_of_rda_files_to_update.log 2>&1 &
###

INFILE <- "rda_files"
OUTFILE <- "rda_files_to_update"

library(BiocInstaller)
available_pkgs <- rownames(available.packages(contrib.url(biocinstallRepos())))

### Install package only if it's not already installed.
installPkg <- function(package)
{
    installed_pkgs <- rownames(installed.packages())
    if (package %in% installed_pkgs)
        return()
    suppressWarnings(suppressPackageStartupMessages(
        library(BiocInstaller, quietly=TRUE)
    ))
    biocLite(package)
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
            cat("  - checking object ", objname, " ... ", sep="")
            obj <- get(objname, envir=envir)
            class_pkg <- attr(class(obj), "package")
            if (is.null(class_pkg) || !(class_pkg %in% available_pkgs)) {
                cat("OK\n")
                next
            }

            ## In order to avoid the "maximal number of DLLs reached..."
            ## infamous error, we use the following trick:
            ##   1) We install 'class_pkg' (if it's not already
            ##      installed) but we do NOT load it.
            ##   2) Then we use a child process to check whether 'obj' is a
            ##      SummarizedExperiment or RangedSummarizedExperiment object.
            ## Unfortunately, this is *much* slower than loading 'class_pkg'
            ## and checking the class of 'obj' in the main process.
            ## Also trying to use detach() didn't prevent to reach the
            ## maximal number of DLLs when we tested this script on
            ## https://hedgehog.fhcrc.org/bioc-data/trunk/experiment/data_store
            ## but maybe I didn't try hard enough (I was only detaching
            ## 'class_pkg' and not all the packages that would get attached
            ## or loaded when doing library(class_pkg)).

            installPkg(class_pkg)

            outline <- paste(rda_path, objname, class(obj), class_pkg, sep="\t")
            Rscript <- c(
                sprintf("suppressWarnings(suppressPackageStartupMessages(library(%s)))", class_pkg),
                sprintf("if (extends('%s', 'SummarizedExperiment')",
                        class(obj)),
                sprintf(" || extends('%s', 'RangedSummarizedExperiment')) {",
                        class(obj)),
                sprintf("    cat('%s\n', file='%s', append=TRUE)",
                        outline, outfile),
                        "    quit(status=0)",
                        "}",
                        "quit(status=1)"
            )
            command <- file.path(R.home("bin"), "R")
            args <- c("--vanilla", "--slave")
            status <- system2(command, args=args, input=Rscript)
            cat(sprintf("OK (status=%d)\n", status))
        }
    }
}

rda_files <- read.table(INFILE, stringsAsFactors=FALSE)[[1L]]
prepareListOfRdaFilesToUpdate(rda_files, outfile=OUTFILE)

