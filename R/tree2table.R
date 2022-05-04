#' @title Convert a NEWICK or NEXUS tree to a table.
#' @description This function converts a tree in NEWICK or NEXUS format to a table (.csv).
#' @param tree.name a file name of the branch length tree.
#' @param type specify the format of branch length tree. The default is NEWICK.
#' @param time logical. If FALSE (default), the input tree is a branch length tree (i.e., not a timetree) and the branch lengths are outputted. If TRUE, the input tree is a timetree and the node times and branch times are outputted.
#' @param filename a file name specified for the output files.
#' @return A table of node label, node ID, descendant nodes, branch lengths, and node times (<filename>.csv).
#' @examples tree2table(tree.name = "example.nwk", type= "NEWICK", filename = "example", time = FALSE)
#' @author Qiqing Tao (qiqing.tao@temple.edu) and Sudhir Kumar
#' @import ape
#' @export

tree2table <- function(tree.name = "", type=c("NEWICK", "NEXUS"), time = FALSE, filename = ""){
  ################# check required packages ##############
  if (!library("ape",logical.return = TRUE)){
    stop("'ape' package not found, please install it to run rrf.times.")
  }

  ## check whether branch length tree is binary
  if (type == "NEXUS"){
    t = ape::read.nexus(tree.name)
  }else{
    t <- ape::read.tree(tree.name)
  }

  if(ape::is.binary(t) == FALSE){
    stop("Only binary trees are allowed. Please remove polytomies.")
  }

  ## convert tree to table ##
  ntips <- ape::Ntip(t)
  if (time == TRUE){
    out.mat <- matrix("-", ncol=7, nrow = ntips+t$Nnode)
    colnames(out.mat) <- c("NodeLabel", "NodeId", "Des1", "Des2", "Brlen1", "Brlen2", "Time")
    out.mat[1:ntips, "NodeLabel"] <- t$tip.label
    out.mat[, "NodeId"] <- c(1:(ntips+t$Nnode))
    out.mat[ntips+c(1:t$Nnode), "Time"] <- ape::branching.times(t)
    for(i in 1:t$Nnode){
      des1 <- which(t$edge[,1]==(ape::Ntip(t)+i))[1]
      des2 <- which(t$edge[,1]==(ape::Ntip(t)+i))[2]
      out.mat[ape::Ntip(t)+i, "Des1"] <- t$edge[des1,2]
      out.mat[ape::Ntip(t)+i, "Des2"] <- t$edge[des2,2]
      out.mat[ape::Ntip(t)+i, "Brlen1"] <- t$edge.length[des1]
      out.mat[ape::Ntip(t)+i, "Brlen2"] <- t$edge.length[des2]
    }
  }else{
    out.mat <- matrix("-", ncol=6, nrow = ntips+t$Nnode)
    colnames(out.mat) <- c("NodeLabel", "NodeId", "Des1", "Des2", "Brlen1", "Brlen2", "Time")
    out.mat[1:ntips, "NodeLabel"] <- t$tip.label
    out.mat[, "NodeId"] <- c(1:(ntips+t$Nnode))
    for(i in 1:t$Nnode){
      des1 <- which(t$edge[,1]==(ape::Ntip(t)+i))[1]
      des2 <- which(t$edge[,1]==(ape::Ntip(t)+i))[2]
      out.mat[ape::Ntip(t)+i, "Des1"] <- t$edge[des1,2]
      out.mat[ape::Ntip(t)+i, "Des2"] <- t$edge[des2,2]
      out.mat[ape::Ntip(t)+i, "Brlen1"] <- t$edge.length[des1]
      out.mat[ape::Ntip(t)+i, "Brlen2"] <- t$edge.length[des2]
    }
  }

  write.csv(data.frame(out.mat), file = paste0(filename, ".csv"), row.names = FALSE)

}
