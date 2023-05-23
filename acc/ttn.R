require(ape, quietly=T)

read.ttn <- function(treefile)
{
    treestr <- readLines(treefile, n=1)
    tree <- read.tree(text=treestr)
    ntips <- Ntip(tree)
    ntipsnodes <- Nnode(tree) + ntips

    all.states <- read.table(treefile, skip=1)
    tip.states <- all.states[1:ntips, 2]
    names(tip.states) <- all.states[1:ntips, 1]
    tip.states <- tip.states[tree$tip.label]

    return(list(tree = tree, states = tip.states))
}
