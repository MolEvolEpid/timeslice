#!/usr/bin/env Rscript

library(ape)
source("tree-stats.R")
source("ttn.R")

get_rules_root <- function(ttnfile)
{
    ttn <- read.ttn(ttnfile)
    tr <- ttn$tree

    tr$tip.label[which(tr$tip.label %in% names(which(ttn$states == 0)))] <- "D"
    tr$tip.label[which(tr$tip.label %in% names(which(ttn$states == 1)))] <- "R"

    dropme <- which(table(tr$edge[,1]) == 1)
    if (length(dropme == 1))
    {
        tr <- diversitree:::drop.tip.fixed(tr, dropme)
    } else if (length(dropme > 1))
    {
        message("HELP!")
    }

    rootstate <- num.intro(tr, "D", "R")$anc
    return(rootstate)
}

# When this file is run as a stand-alone script:
if (sys.nframe() == 0)
{
    wd <- commandArgs(trailingOnly=TRUE)[1]
    setwd(wd)
    ttnfiles <- Sys.glob("*.ttn")

    ans <- sapply(ttnfiles, get_rules_root)
    names(ans) <- basename(ttnfiles)

    write.table(ans, file="rules.csv", quote=F, col.names=F, sep=",")
    table(ans)
}
