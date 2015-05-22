### =========================================================================
### prepare_list_of_rda_files_to_update.R
### -------------------------------------------------------------------------
###
### This script performs STEP 3 of the "Find and update objects" procedure
### described in the 'Find_and_update_objects.txt' file located in the same
### folder.
###
### Before you run this script, make sure you performed STEPS 1 & 2.
### See Find_and_update_objects.txt for more information.
###
### Then to run STEP 3 in "batch" mode:
###
###   cd <dir/you/want/to/search>  # the 'rda_objects' file should be here
###   R CMD BATCH prepare_list_of_rda_files_to_update.R \
###              >prepare_list_of_rda_files_to_update.log 2>&1 &
###

INFILE <- "rda_objects"
OUTFILE <- "rda_files_to_update"

library(BiocInstaller)
library(SummarizedExperiment)

prepareListOfRdaFilesToUpdate <- function(rda_objects, outfile="")
{
    rda_objects2 <- unique(rda_objects[ , c("objclass", "objclass_pkg")])
    objclass2 <- rda_objects2[ , "objclass"]
    objclass_pkg2 <- rda_objects2[ , "objclass_pkg"]
    idx <- which(duplicated(objclass2))
    if (length(idx) != 0L) {
        msg <- c("the following classes are defined in more than 1 package: ",
                 paste0(unique(objclass2[idx]), collapse=", "))
        warning(msg)
    }
    pkg2class <- split(objclass2, objclass_pkg2)
    pkg2class[c(".", ".GlobalEnv")] <- NULL
    pkgs <- names(pkg2class)

    ## Install missing packages.
    installed_pkgs <- rownames(installed.packages())
    missing_pkgs <- setdiff(pkgs, installed_pkgs)
    if (length(missing_pkgs) != 0L) {
        biocLite(missing_pkgs)
        installed_pkgs <- rownames(installed.packages())
        missing_pkgs <- setdiff(pkgs, installed_pkgs)
        if (length(missing_pkgs) != 0L) {
            ## Some packages could not be installed.
            pkgs <- intersect(pkgs, installed_pkgs)
            pkg2class <- pkg2class[pkgs]
        }
    }

    ## Load packages where all object classes are defined.
    ## Unfortunately, this sometimes hits the maximal number of DLLs that can
    ## be loaded ("maximal number of DLLs reached..." error).
    #for (i in seq_along(pkgs)) {
    #    pkg <- pkgs[[i]]
    #    cat("[", i , "/", length(pkgs), "] Loading ", pkg, " ... ",
    #        sep="")
    #    suppressWarnings(suppressPackageStartupMessages(
    #        library(pkg, character.only=TRUE, quietly=TRUE)
    #    ))
    #    cat("OK\n")
    #}

    ## We use subprocesses to work around the "maximal number of DLLs
    ## reached..." error.
    for (i in seq_along(pkgs)) {
        pkg <- pkgs[[i]]
        cat("[", i , "/", length(pkgs), "] Check classes defined in package ",
            pkg, " ... ", sep="")
        classes <- paste0("\"", pkg2class[[pkg]], "\"")
        classes <- paste0("c(", paste(classes, collapse=", "), ")")
        class_summary_file <- paste0(pkg, "_class_summary")
        input <- c(
                "suppressWarnings(suppressPackageStartupMessages(",
        sprintf("    library(%s)", "SummarizedExperiment"),
                "))",
                "suppressWarnings(suppressPackageStartupMessages(",
        sprintf("    library(%s)", pkg),
                "))",
        sprintf("ok <- sapply(%s, function(class) {", classes),
                "    extends(class, \"SummarizedExperiment\") ||",
                "    extends(class, \"RangedSummarizedExperiment\")",
                "})",
                "class_summary <- data.frame(class=names(ok), ok=unname(ok))",
        sprintf("write.table(class_summary, file=\"%s\",", class_summary_file),
                "            sep=\"\t\")"
        )
        command <- file.path(R.home("bin"), "R")
        args <- c("--vanilla", "--slave")
        system2(command, args=args, input=input)
        cat("OK\n")
    }

    ## Collect the class summary files.
    class_summaries <- lapply(pkgs, function(pkg) {
        class_summary_file <- paste0(pkg, "_class_summary")
        read.table(class_summary_file, stringsAsFactors=FALSE)
    })
    class_summary <- do.call(rbind, class_summaries)

    ## Sanity check.
    class1 <- class_summary[ , "class"]
    class2 <- unlist(pkg2class, use.names=FALSE)
    stopifnot(identical(class1, class2))

    ## Write output to file.
    class2ok <- setNames(class_summary[ , "ok"], class1)
    objclass <- rda_objects[ , "objclass"]
    ok <- class2ok[objclass]
    ok[is.na(ok)] <- FALSE
    rda_files_to_update <- rda_objects[ok, , drop=FALSE]
    rda_files_to_update <- do.call(
        paste,
        c(as.list(rda_files_to_update), list(sep="\t"))
    )
    writeLines(rda_files_to_update, con=outfile)
}

rda_objects <- read.table(INFILE, stringsAsFactors=FALSE)
colnames(rda_objects) <- c("rda_path", "objname", "objclass", "objclass_pkg")
prepareListOfRdaFilesToUpdate(rda_objects, outfile=OUTFILE)

