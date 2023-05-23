#!/usr/bin/env Rscript

require(stringr)
source("ttn.R")

get_topology <- function(treefile)
{
    capture.output(ttn <- read.ttn(treefile), file="/dev/null")

    ans <- "PP"
    tips <- names(ttn$states[ttn$states == 0])
    if (is.monophyletic(ttn$tree, tips))
        stringr::str_sub(ans, 1, 1) <- "M"
    tips <- names(ttn$states[ttn$states == 1])
    if (is.monophyletic(ttn$tree, tips))
        stringr::str_sub(ans, 2, 2) <- "M"

    return(ans)
}

# When this file is run as a stand-alone script:
if (sys.nframe() == 0)
{
    wd <- commandArgs(trailingOnly=TRUE)[1]
    setwd(wd)
    ttnfiles <- Sys.glob("*.ttn")

    ans <- sapply(ttnfiles, get_topology)

    write.table(ans, file="topology.csv", quote=F, col.names=F, sep=",")
    table(ans)
}
