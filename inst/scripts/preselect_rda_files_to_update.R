### =========================================================================
### preselect_rda_files_to_update.R
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
###   R CMD BATCH preselect_rda_files_to_update.R \
###              >preselect_rda_files_to_update.log 2>&1 &
###

INFILE <- "rda_objects"
OUTFILE <- "rda_files_to_update"

library(BiocInstaller)
library(SummarizedExperiment)

if (FALSE) {
### Unfortunately, loading all the required packages in the main process will
### sometimes hit the maximal number of DLLs that can be loaded ("maximal
### number of DLLs reached..." infamous error).
checkClasses <- function(classes, package)
{
    suppressWarnings(suppressPackageStartupMessages(
        library(package, character.only=TRUE, quietly=TRUE)
    ))
    sapply(classes, function(class) {
           extends(class, "SummarizedExperiment") ||
           extends(class, "RangedSummarizedExperiment")
    })
}
}

### We check the classes in a subprocess to work around the "maximal number
### of DLLs reached..." infamous error.
checkClasses <- function(classes, package)
{
    classes_in1string <- paste0("\"", classes, "\"")
    classes_in1string <- paste0("c(",
                                paste(classes_in1string, collapse=", "),
                                ")")
    outfile <- file.path(tempdir(), paste0(package, "_class_summary"))
    input <- c("suppressWarnings(suppressPackageStartupMessages(",
       sprintf("    library(%s)", "SummarizedExperiment"),
               "))",
               "suppressWarnings(suppressPackageStartupMessages(",
       sprintf("    library(%s)", package),
               "))",
       sprintf("classes <- %s", classes_in1string),
               "ok <- sapply(classes, function(class) {",
               "    extends(class, \"SummarizedExperiment\") ||",
               "    extends(class, \"RangedSummarizedExperiment\")",
               "})",
               "class_summary <- data.frame(class=classes, ok=unname(ok))",
       sprintf("write.table(class_summary, file=\"%s\", sep=\"\t\")", outfile)
    )
    command <- file.path(R.home("bin"), "R")
    args <- c("--vanilla", "--slave")
    system2(command, args=args, input=input)
    class_summary <- read.table(outfile, stringsAsFactors=FALSE)
    file.remove(outfile)
    stopifnot(identical(class_summary[ , "class"], classes))  # sanity check
    setNames(class_summary[ , "ok"], classes)
}

preselectRdaFilesToUpdate <- function(rda_objects, outfile="")
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

    class2ok <- unlist(
        lapply(seq_along(pkgs),
            function(i) {
                pkg <- pkgs[[i]]
                cat("[", i , "/", length(pkgs), "] Check classes defined ",
                    "in package ", pkg, " ... ", sep="")
                ans <- checkClasses(pkg2class[[pkg]], pkg)
                cat("OK\n")
                ans
            }),
        use.names=FALSE
    )

    ## Write output to file.
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
preselectRdaFilesToUpdate(rda_objects, outfile=OUTFILE)

