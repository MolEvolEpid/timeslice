#!/usr/bin/env Rscript

require(phytools)
require(phangorn)

#--------------------------------------------------
# Helper functions
#--------------------------------------------------

# For reading in trees with states.
source("ttn.R")

# To a phylo object, add node.time element that is the distance from each tip or
# node to the root.
add_node_times <- function(phy)
{
    # order = nodes out from root, then tips
    i_node <- c((Ntip(phy)+1) : (Ntip(phy)+Nnode(phy)), seq(Ntip(phy)))
    phy$node.time <- rep(NA, length(i_node))

    # at the root
    phy$node.time[i_node[1]] <- 0

    # all other nodes and tips
    for (n in i_node[-1])
        phy$node.time[n] <- phy$edge.length[which(phy$edge[,2] == n)] + phy$node.time[phytools::getParent(phy, n)]

    return(phy)
}

#--------------------------------------------------
# Main function
#--------------------------------------------------

# The first value is the earliest time transmission could logically have occurred.
#   It's the most recent node that has daughters in both states.
#   (So actually, the earliest time is immediately after this time.)
# The second value is the time of the latest sample.

# NOTE: The end of the possible window is the first time by which there exists
# a sample from each person.  So if someone is sampled at more than one time,
# the possible window could end before the final sample date on the tree.
# But the final sample date is what's reported.

possible_transmission_window <- function(treefile)
{
    # Read in the tree and states
    capture.output(ttn <- read.ttn(treefile), file="/dev/null")
    ttn$tree <- add_node_times(ttn$tree)
    phy <- ttn$tree

    # Get the order in which to work through the nodes (shallow before deep)
    # early elements are deep nodes, etc:
    nodes_ordered <- order(phy$node.time)
    # drop the ones that are tips:
    nodes_ordered <- nodes_ordered[nodes_ordered > Ntip(phy)]
    # put the latest nodes first:
    nodes_ordered <- rev(nodes_ordered)

    # Work through the nodes in order
    for (node in nodes_ordered)
    {
        i_tips <- phangorn::Descendants(phy, node, "tips")[[1]]
        tip_states <- ttn$states[i_tips]
        if (all(c(0, 1) %in% tip_states))
           break
    }
    early <- phy$node.time[node]
    late <- max(phy$node.time)

    return(c(early, late))
}

# When this file is run as a stand-alone script:
if (sys.nframe() == 0)
{
    input <- commandArgs(trailingOnly=TRUE)[1]

    if (file.info(input)$isdir)
    {
        # Operate on a whole directory of files
        setwd(input)
        ttnfiles <- Sys.glob("*.ttn")

        ans <- t(sapply(ttnfiles, possible_transmission_window))

        write.table(ans, file="window.csv", quote=F, col.names=F, sep=",")

    } else {
        # Operate on a single file
        ans <- possible_transmission_window(input)
        cat(ans, "\n")
    }
}
