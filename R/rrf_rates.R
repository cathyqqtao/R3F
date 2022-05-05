#' @title Estimate relative lineage rates.
#' @description This function estimates relative lineage rates from a branch length tree using the relative rate framework (RRF).
#' @param tree.name a file name of the branch length tree.
#' @param type specify the format of branch length tree. The default is NEWICK.
#' @param outgroup a character list containing names of tips belonging to the rooting outgroup, which will be removed in the calculation. If outgroup = "" (the default), the input tree must be rooted and no tips will be removed.
#' @param filename a file name specified for the output file.
#' @return A table of relative lineage rates (<filename>_RRF_rates.csv).
#' @return A tree with relative rates on branches (i.e., branch lengths are rates) (<filename>_RRF_rates.nwk).
#' @examples rrf_rates(tree.name = "example.nwk", type= "NEWICK", outgroup = c("Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus"), filename = "example")
#' @author Qiqing Tao (qiqing.tao@temple.edu) and Sudhir Kumar
#' @references K. Tamura et al. Mol. Biol. Evol. (2018) 35:1770-1782. doi:10.1093/molbev/msy044.
#' @import ape
#' @importFrom phangorn Descendants
#' @export
#'
rrf_rates <- function(tree.name = "", type=c("NEWICK", "NEXUS"), outgroup = "", filename = ""){

  ################# check required packages ##############
  if (!library("ape",logical.return = TRUE)){
    stop("'ape' package not found, please install it to run rrf.rates.")
  }
  if (!library("phangorn",logical.return = TRUE)){
    stop("'phangorn' package not found, please install it to run rrf.rates.")
  }

  ################# check branch length tree and outgroup #########
  if (type == "NEXUS"){
    t <- ape::read.nexus(tree.name)
  }else{
    t <- ape::read.tree(tree.name)
  }

  ## check outgroups
  suppressWarnings(if (outgroup != ""){
    for (i in 1:length(outgroup)){
      if(is.na(match(outgroup[i], t$tip.label)) == TRUE){
        stop(paste("Outgroup \"",outgroup[i], "\" is not found. Please check.",sep=''))
      }
    }
  })

  ## check whether branch length tree is binary
  if(ape::is.binary(t) == FALSE){
    stop("Only binary trees are allowed. Please remove polytomies.")
  }


  ########################################################################
  ###################### calculate relative lineage rates ################
  ########################################################################

  ## check whether the tree is rooted
  suppressWarnings(if (outgroup != ""){
    t <- ape::root(t, outgroup, resolve.root = TRUE)
    t <- ape::drop.tip(t, outgroup)
  }else{
    if (is.rooted(t) == FALSE){
      stop("Please provide a rooted tree or the names of tips in the rooting outgroup.")
    }
  })

  ## check whether there are too many zero-length branches
  brlen <- t$edge.length
  zero.brlen = length(which(brlen==0))/length(brlen)
  if (zero.brlen >= 0.1){
    print(paste0("Warning: there are too many zero-length braches (>", floor(zero.brlen*100), "%) in the tree, rates calculation may not be reliable on those branches."))
  }

  #### get raw relative rates and node times ####

  RRF.mat <- matrix(data=NA, nrow = t$Nnode, ncol = 13)
  colnames(RRF.mat) <- c('NodeId', 'Des1', 'Des2', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'r5', 'r6', 'r5.adjust', 'r6.adjust')

  tips.num <- length(t$tip.label)
  RRF.mat[, 1] <- seq(tips.num+1, tips.num+t$Nnode, 1)

  for (nd in (tips.num+t$Nnode):(tips.num+1)){  # from the shallowest internal node to the root
    row.id <- match(nd, RRF.mat[, 1])

    des1 <- phangorn::Descendants(t, nd, 'children')[1]
    des2 <- phangorn::Descendants(t, nd, 'children')[2]
    RRF.mat[row.id, 2] <- as.integer(des1)
    RRF.mat[row.id, 3] <- as.integer(des2)


    if (des1 <= tips.num && des2 <= tips.num){ ## 2-clade case = des1 and 2 are tips
      l1 <- 0
      l2 <- 0
      l3 <- 0
      l4 <- 0
      l5 <- brlen[match(des1, t$edge[, 2])]
      l6 <- brlen[match(des2, t$edge[, 2])]

      RRF.mat[row.id, seq(4,9,1)] <- c(l1,l2,l3,l4,l5,l6)

    }


    if (des1 > tips.num && des2 > tips.num){ ## 4-clade case
      des1.row.id <- match(des1, RRF.mat[, 1])
      des2.row.id <- match(des2, RRF.mat[, 1])

      l1 <- sqrt(RRF.mat[des1.row.id, 4]*RRF.mat[des1.row.id, 5])+RRF.mat[des1.row.id, 8]
      l2 <- sqrt(RRF.mat[des1.row.id, 6]*RRF.mat[des1.row.id, 7])+RRF.mat[des1.row.id, 9]
      l3 <- sqrt(RRF.mat[des2.row.id, 4]*RRF.mat[des2.row.id, 5])+RRF.mat[des2.row.id, 8]
      l4 <- sqrt(RRF.mat[des2.row.id, 6]*RRF.mat[des2.row.id, 7])+RRF.mat[des2.row.id, 9]
      l5 <- brlen[match(des1, t$edge[,2])]
      l6 <- brlen[match(des2, t$edge[,2])]

      RRF.mat[row.id, seq(4,9,1)] <- c(l1,l2,l3,l4,l5,l6)

      if (l1 == 0){ ## to avoid INF problem in rate calculation
        l1 <- 10e-20
      }
      if (l2 == 0){
        l2 <- 10e-20
      }
      if (l3 == 0){
        l3 <- 10e-20
      }
      if (l4 == 0){
        l4 <- 10e-20
      }
      if (sqrt(l3*l4)+l6 == 0){ ## to avoid INF problem in rate calculation
        l6 <- 1e-20
      }
      if (sqrt(l1*l2)+l5 == 0){
        l5 <- 1e-20
      }

      r5 <- sqrt(sqrt(l1*l2)+l5)/sqrt(sqrt(l3*l4)+l6)
      r6 <- sqrt(sqrt(l3*l4)+l6)/sqrt(sqrt(l1*l2)+l5)

      RRF.mat[row.id, c(10,11)] <- c(r5,r6)

      des1.des1 <- phangorn::Descendants(t, des1, "children")[1]
      des1.des2 <- phangorn::Descendants(t, des1, "children")[2]
      des2.des1 <- phangorn::Descendants(t, des2, "children")[1]
      des2.des2 <- phangorn::Descendants(t, des2, "children")[2]

      r1 <- sqrt(l1)*sqrt(sqrt(l1*l2)+l5)/(sqrt(l2)*sqrt(sqrt(l3*l4)+l6))
      r2 <- sqrt(l2)*sqrt(sqrt(l1*l2)+l5)/(sqrt(l1)*sqrt(sqrt(l3*l4)+l6))
      r3 <- sqrt(l3)*sqrt(sqrt(l3*l4)+l6)/(sqrt(l4)*sqrt(sqrt(l1*l2)+l5))
      r4 <- sqrt(l4)*sqrt(sqrt(l3*l4)+l6)/(sqrt(l3)*sqrt(sqrt(l1*l2)+l5))

      if (des1.des1 <= tips.num){ ## If any 4 clade contains only 1 tip, we need to use the calculated r1,r2,r3,r4 as the tip rates
        RRF.mat[des1.row.id, 10] = r1
      }
      if (des1.des2 <= tips.num){
        RRF.mat[des1.row.id, 11] = r2
      }
      if (des2.des1 <= tips.num){
        RRF.mat[des2.row.id, 10] = r3
      }
      if (des2.des2 <= tips.num){
        RRF.mat[des2.row.id, 11] = r4
      }

    }


    if (des1 > tips.num && des2 <= tips.num){ ## 3-clade, des2 is the tip
      des1.row.id <- match(des1, RRF.mat[, 1])

      l1 <- sqrt(RRF.mat[des1.row.id, 4]*RRF.mat[des1.row.id, 5])+RRF.mat[des1.row.id, 8]
      l2 <- sqrt(RRF.mat[des1.row.id, 6]*RRF.mat[des1.row.id, 7])+RRF.mat[des1.row.id, 9]
      l3 <- 0
      l4 <- 0
      l5 <- brlen[match(des1, t$edge[, 2])]
      l6 <- brlen[match(des2, t$edge[, 2])]

      RRF.mat[row.id, seq(4,9,1)] <- c(l1,l2,l3,l4,l5,l6)

      if (l1 == 0){ ## to avoid INF problem in rate calculation
        l1 <- 10e-20
      }
      if (l2 == 0){
        l2 <- 10e-20
      }
      if (l6 == 0){ ## to avoid INF problem in rate calculation
        l6 <- 10e-20
      }

      r5 <- sqrt(sqrt(l1*l2)+l5)/sqrt(l6)
      r6 <- sqrt(l6)/sqrt(sqrt(l1*l2)+l5)

      RRF.mat[row.id, c(10,11)] <- c(r5,r6)

      des1.des1 <- phangorn::Descendants(t, des1, "children")[1]
      des1.des2 <- phangorn::Descendants(t, des1, "children")[2]

      r1 <- sqrt(l1)*sqrt(sqrt(l1*l2)+l5)/sqrt(l2*l6)
      r2 <- sqrt(l2)*sqrt(sqrt(l1*l2)+l5)/sqrt(l1*l6)

      if (des1.des1 <= tips.num){
        RRF.mat[des1.row.id, 10] <- r1
      }
      if (des1.des2 <= tips.num){
        RRF.mat[des1.row.id, 11] = r2
      }

    }


    if (des1 <= tips.num && des2 > tips.num){ ## 3-clade, des1 is the tip
      des2.row.id <-match(des2, RRF.mat[, 1])

      l1 <- 0
      l2 <- 0
      l3 <- sqrt(RRF.mat[des2.row.id, 4]*RRF.mat[des2.row.id, 5])+RRF.mat[des2.row.id, 8]
      l4 <- sqrt(RRF.mat[des2.row.id, 6]*RRF.mat[des2.row.id, 7])+RRF.mat[des2.row.id, 9]
      l5 <- brlen[match(des1, t$edge[, 2])]
      l6 <- brlen[match(des2, t$edge[, 2])]

      RRF.mat[row.id, seq(4,9,1)] <- c(l1,l2,l3,l4,l5,l6)

      if (l5 == 0){ ## to avoid INF problem in rate calculation
        l5 <- 10e-20
      }
      if (l3 == 0){ ## to avoid INF problem in rate calculation
        l3 <- 10e-20
      }
      if (l4 == 0){
        l4 <- 10e-20
      }

      r5 <- sqrt(l5)/sqrt(sqrt(l3*l4)+l6)
      r6 <- sqrt(sqrt(l3*l4)+l6)/sqrt(l5)

      RRF.mat[row.id, c(10,11)] <- c(r5,r6)

      des2.des1 <- phangorn::Descendants(t, des2, "children")[1]
      des2.des2 <- phangorn::Descendants(t, des2, "children")[2]

      r3 <- sqrt(l3)*sqrt(sqrt(l3*l4)+l6)/sqrt(l4*l5)
      r4 <- sqrt(l4)*sqrt(sqrt(l3*l4)+l6)/sqrt(l3*l5)

      if (des2.des1 <= tips.num){
        RRF.mat[des2.row.id, 10] <- r3
      }
      if (des2.des2 <= tips.num){
        RRF.mat[des2.row.id, 11] <- r4
      }

    }
  }

  # RRF.mat[is.na(RRF.mat[,10]), 10] <- 0
  # RRF.mat[is.infinite(RRF.mat[,10]), 10] <- 0
  # RRF.mat[is.na(RRF.mat[,11]), 11] <- 0
  # RRF.mat[is.infinite(RRF.mat[,11]), 11] <- 0

  #### adjust relative rates by multiplying the ancestral rate ####
  for (nd in (tips.num+1):(tips.num+t$Nnode)){ ## from root to shallowest internal nodes

    row.id <- match(nd, RRF.mat[, 1])
    anc.row <- c(match(nd, RRF.mat[,2]), match(nd, RRF.mat[,3]))

    if (is.na(anc.row[1]) == TRUE && is.na(anc.row[2]) == TRUE){ ## ancestor is root
      r.anc <- 1
      RRF.mat[row.id, c(12,13)] <- RRF.mat[row.id, c(10,11)]*r.anc
    }else if (is.na(anc.row[1]) == FALSE && is.na(anc.row[2]) == TRUE){
      r.anc <- RRF.mat[anc.row[1], 12]
      RRF.mat[row.id, c(12,13)] <- RRF.mat[row.id, c(10,11)]*r.anc

      if(RRF.mat[row.id, 2] <= tips.num){ ## if one offspring is tip, we need grandparent rate
        grandpa.row <- c(match(RRF.mat[anc.row[1], 1], RRF.mat[, 2]), match(RRF.mat[anc.row[1], 1], RRF.mat[, 3]))
        if (is.na(grandpa.row[1]) == TRUE && is.na(grandpa.row[2]) == TRUE){
          r.anc <- 1
        }else{
          r.grandpa <- c(RRF.mat[grandpa.row[1], 12], RRF.mat[grandpa.row[2], 13])
          r.anc <-  r.grandpa[!is.na(r.grandpa)]
          RRF.mat[row.id, 12] <- RRF.mat[row.id, 10]*r.anc
        }
      }

      if(RRF.mat[row.id, 3] <= tips.num){ ## if one offspring is tip, we need grandparent rate
        grandpa.row <- c(match(RRF.mat[anc.row[1], 1], RRF.mat[, 2]), match(RRF.mat[anc.row[1], 1], RRF.mat[, 3]))
        if (is.na(grandpa.row[1]) == TRUE && is.na(grandpa.row[2]) == TRUE){
          r.anc <- 1
        }else{
          r.grandpa <- c(RRF.mat[grandpa.row[1], 12], RRF.mat[grandpa.row[2], 13])
          r.anc <-  r.grandpa[!is.na(r.grandpa)]
          RRF.mat[row.id, 13] <- RRF.mat[row.id, 11]*r.anc
        }
      }
    }else if (is.na(anc.row[1]) == TRUE && is.na(anc.row[2]) == FALSE){
      r.anc <- RRF.mat[anc.row[2], 13]
      RRF.mat[row.id, c(12,13)] <- RRF.mat[row.id, c(10,11)]*r.anc

      if(RRF.mat[row.id, 2] <= tips.num){ ## if one offspring is tip, we need grandparent rate
        grandpa.row <- c(match(RRF.mat[anc.row[2], 1], RRF.mat[, 2]), match(RRF.mat[anc.row[2], 1], RRF.mat[, 3]))
        if (is.na(grandpa.row[1]) == TRUE && is.na(grandpa.row[2]) == TRUE){
          r.anc <- 1
        }else{
          r.grandpa <- c(RRF.mat[grandpa.row[1], 12], RRF.mat[grandpa.row[2], 13])
          r.anc <-  r.grandpa[!is.na(r.grandpa)]
          RRF.mat[row.id, 12] <- RRF.mat[row.id, 10]*r.anc
        }
      }

      if(RRF.mat[row.id, 3] <= tips.num){ ## if one offspring is tip, we need grandparent rate
        grandpa.row <- c(match(RRF.mat[anc.row[2], 1], RRF.mat[, 2]), match(RRF.mat[anc.row[2], 1], RRF.mat[, 3]))
        if (is.na(grandpa.row[1]) == TRUE && is.na(grandpa.row[2]) == TRUE){
          r.anc <- 1
        }else{
          r.grandpa <- c(RRF.mat[grandpa.row[1], 12], RRF.mat[grandpa.row[2], 13])
          r.anc <-  r.grandpa[!is.na(r.grandpa)]
          RRF.mat[row.id, 13] <- RRF.mat[row.id, 11]*r.anc
        }
      }
    }
  }

  ############ END OF CALCULATION ##############

  ####### generate output table ########
  out.mat <- matrix(NA, ncol=5, nrow = t$Nnode+tips.num)
  colnames(out.mat) <- c("NodeLabel", "NodeId", "Des1", "Des2", "Rate")
  out.mat[, "NodeLabel"] <- c(t$tip.label, rep("-", t$Nnode))
  out.mat[, "NodeId"] <- c(1:(t$Nnode+tips.num))
  out.mat[, "Des1"] <- c(rep("-", tips.num), RRF.mat[,"Des1"])
  out.mat[, "Des2"] <- c(rep("-", tips.num), RRF.mat[,"Des2"])
  out.mat[RRF.mat[,"Des1"], "Rate"] <- RRF.mat[, "r5.adjust"]
  out.mat[RRF.mat[,"Des2"], "Rate"] <- RRF.mat[, "r6.adjust"]
  out.mat[tips.num+1, "Rate"] <- sprintf("%1.3f",1) ## root rate

  write.csv(data.frame(out.mat), file = paste0(filename, "_RRF_rates.csv"), row.names = FALSE)

  ####### generate output tree (nexus) with rates ########
  out.tree <- t
  out.tree$edge.length <- as.numeric(out.mat[out.tree$edge[,2], "Rate"])

  ape::write.tree(out.tree, file = paste0(filename, "_RRF_rates.nwk"))
}
