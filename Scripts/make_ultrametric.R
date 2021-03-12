#---------------------------------------------------------------------------------------------------
#MAKE ULTRAMETRIC TREES
#---------------------------------------------------------------------------------------------------

make.ultrametric.tree <- function(tree) {
    root.number <- length(tree$tip.label) + 1
    ultra.tree <- length.normalise(tree, tree, root.number, 1)
    return(ultra.tree)
}

length.normalise <- function(orig.tree, new.tree, curr.node, remaining.stick) {
    curr.node.children <- unlist(Descendants(orig.tree, curr.node, "children"))

    for (j in curr.node.children) {
      index <- which(orig.tree$edge[,2] == j & orig.tree$edge[,1] == curr.node, arr.ind = TRUE)

      if (j %in% orig.tree$tip.label) {
        new.tree$edge.length[index] <- remaining.stick
      } else {
        curr.node.tips <- unlist(Descendants(orig.tree, j, "tips"))
        curr.dist <- find.distance(orig.tree, curr.node, j)
        if (curr.dist == 0) {curr.dist <- 0.01} # So that no edge lengths are zero
        desc.lengths <- sapply(curr.node.tips, FUN = find.distance, tree = orig.tree, from = curr.node)
        if (mean(desc.lengths) != 0) {
          new.tree$edge.length[index] <- remaining.stick * curr.dist / mean(desc.lengths)
        } else {
          new.tree$edge.length[index] <- remaining.stick
        }
        shorter.stick <- remaining.stick - new.tree$edge.length[index]

        # Call as recursive function
        new.tree <- length.normalise(orig.tree, new.tree, j, shorter.stick)
      }
    }
    return(new.tree)
} 
  
find.distance <- function(tree, from, to) {
  path <- nodepath(tree, from, to)
  res <- 0
  for (j in 2:length(path)) {
    index <- which(tree$edge[,2] == path[j] & tree$edge[,1] == path[j-1], arr.ind = TRUE)
    res <- res + tree$edge.length[index]
  }
  return(res)
}  
