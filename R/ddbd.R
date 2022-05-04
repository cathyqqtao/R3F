#' @title Estimate parameter values for Birth-Death speciation prior.
#' @description This function estimates the values of parameters of Birth-Death speciation tree prior for downstream Bayesian dating analysis from a branch length tree.
#' @param tree.name a file name of the branch length tree.
#' @param type specify the format of branch length tree. The default is NEWICK.
#' @param outgroup a character list containing names of tips belonging to the rooting outgroup, which will be removed in the calculation. If outgroup = "" (the default), the input tree must be rooted and no tips will be removed.
#' @param sampling.frac a numeric value representing the sampling fraction. If sampling.frac = 0, sampling fraction will be estimated. Otherwise, it will be fixed to the given value.
#' @param anchor.node an ingeter corresponding to the ID of a node with user-provide node time. The user-provide value is used to convert relative rates to absolute rates. If anchor.node = 0 (the default), rates will not be converted and are still relative rates.
#' @param anchor.time a numeric value specifying to the node time of anchor.node. When anchor.node = 0, anchor.time will be ingnored. See "details".
#' @param measure specify the method for selecting the initial values. The best initial values can be selected by minimzing the sum of squared errors (SSE) or Kullback-Leibler divergence (KL). The default is SSE.
#' @param filename a file name specified for the output file.
#' @param plot logical. If TRUE (default), a histogram of node times and a density curve of the estimated parameters of Birth-Death model are plotted.
#' @details Birth rate and death rate are estimated using the same time unit as provided in "anchor.time". It is recommened to adjust the time unit, so that the maximum value of "anchor time" does not exceed 10.
#' @return Values of parameters in Birth-Death speciation model (<filename>_ddbd.txt). If the sampling.fraction is specified, only birth rate and death rate will be estimated. Otherwise, all three parameters will be estimated.
#' @examples ddbd(tree.name = "example.nwk", type= "NEWICK", outgroup = c("Ornithorhynchus_anatinus", "Zaglossus_bruijni", "Tachyglossus_aculeatus"), anchor.node = 272, anchor.time = 1.85, measure = "SSE", filename = "example", plot = TRUE)
#' @author Qiqing Tao (qiqing.tao@temple.edu) and Sudhir Kumar
#' @references Q. Tao et al. Bioinformatics (2021) 37:i102-i110. doi:10.1093/bioinformatics/btab307
#' @import ape
#' @importFrom phangorn Descendants Ancestors
#' @importFrom stats4 mle
#' @importFrom FNN KL.divergence
#' @export

ddbd <- function(tree.name = "", type=c("NEWICK", "NEXUS"), outgroup = "", sampling.frac = 0, anchor.node = 0, anchor.time = 1, measure = c("SSE","KL"), filename = "", plot = TRUE){

  ################# check required packages ##############
  if (!library("ape",logical.return = TRUE)){
    stop("'ape' package not found, please install it to run ddbd.")
  }
  if (!library("phangorn",logical.return = TRUE)){
    stop("'phangorn' package not found, please install it to run ddbd.")
  }
  if (!library("stats4",logical.return = TRUE)){
    stop("'stats4' package not found, please install it to run ddbd.")
  }
  if (!library("FNN",logical.return = TRUE)){
    stop("'FNN' package not found, please install it to run ddbd.")
  }

  ################# check brach length tree and outgroup #########
  if (type == "NEXUS"){
    t = ape::read.nexus(tree.name)
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
  ######### calculate relative lineage rates and node times ##############
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
    print(paste0("Warning: there are too many zero-length braches (>", floor(zero.brlen*100), "%) in the tree, rates and node times calculation may not be reliable for nodes connecting to those branches."))
  }

  #### get raw relative rates and node times ####

  RRF.mat <- matrix(data=NA, nrow = t$Nnode, ncol = 15)
  colnames(RRF.mat) <- c('NodeId', 'Des1', 'Des2', 'l1', 'l2', 'l3', 'l4', 'l5', 'l6', 'r5', 'r6', 'r5.adjust', 'r6.adjust', "t7", "t7.adjust")

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

      t7 <- sqrt(sqrt(l1*l2)+l5)*sqrt(sqrt(l3*l4)+l6)
      RRF.mat[row.id, "t7"] <- t7

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
        t5 <- sqrt(l1*l2)*sqrt(sqrt(l3*l4)+l6)/sqrt(sqrt(l1*l2)+l5)
        RRF.mat[des1.row.id, "t7"] <- t5
      }
      if (des1.des2 <= tips.num){
        RRF.mat[des1.row.id, 11] = r2
        t5 <- sqrt(l1*l2)*sqrt(sqrt(l3*l4)+l6)/sqrt(sqrt(l1*l2)+l5)
        RRF.mat[des1.row.id, "t7"] <- t5
      }
      if (des2.des1 <= tips.num){
        RRF.mat[des2.row.id, 10] = r3
        t6 <- sqrt(l3*l4)*sqrt(sqrt(l1*l2)+l5)/sqrt(sqrt(l3*l4)+l6)
        RRF.mat[des2.row.id, "t7"] <- t6
      }
      if (des2.des2 <= tips.num){
        RRF.mat[des2.row.id, 11] = r4
        t6 <-sqrt(l3*l4)*sqrt(sqrt(l1*l2)+l5)/sqrt(sqrt(l3*l4)+l6)
        RRF.mat[des2.row.id, "t7"] <- t6
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

      t7 <- sqrt(sqrt(l1*l2)+l5)*sqrt(sqrt(l3*l4)+l6)
      RRF.mat[row.id, "t7"] <- t7

      des1.des1 <- phangorn::Descendants(t, des1, "children")[1]
      des1.des2 <- phangorn::Descendants(t, des1, "children")[2]

      r1 <- sqrt(l1)*sqrt(sqrt(l1*l2)+l5)/sqrt(l2*l6)
      r2 <- sqrt(l2)*sqrt(sqrt(l1*l2)+l5)/sqrt(l1*l6)

      if (des1.des1 <= tips.num){
        RRF.mat[des1.row.id, 10] <- r1
        t5 <- sqrt(l1*l2)*sqrt(sqrt(l3*l4)+l6)/sqrt(sqrt(l1*l2)+l5)
        RRF.mat[des1.row.id, "t7"] <- t5
      }
      if (des1.des2 <= tips.num){
        RRF.mat[des1.row.id, 11] = r2
        t5 <- sqrt(l1*l2)*sqrt(sqrt(l3*l4)+l6)/sqrt(sqrt(l1*l2)+l5)
        RRF.mat[des1.row.id, "t7"] <- t5
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

      t7 <- sqrt(sqrt(l1*l2)+l5)*sqrt(sqrt(l3*l4)+l6)
      RRF.mat[row.id, "t7"] <- t7

      des2.des1 <- phangorn::Descendants(t, des2, "children")[1]
      des2.des2 <- phangorn::Descendants(t, des2, "children")[2]

      r3 <- sqrt(l3)*sqrt(sqrt(l3*l4)+l6)/sqrt(l4*l5)
      r4 <- sqrt(l4)*sqrt(sqrt(l3*l4)+l6)/sqrt(l3*l5)

      if (des2.des1 <= tips.num){
        RRF.mat[des2.row.id, 10] <- r3
        t6 <- sqrt(l3*l4)*sqrt(sqrt(l1*l2)+l5)/sqrt(sqrt(l3*l4)+l6)
        RRF.mat[des2.row.id, "t7"] <- t6
      }
      if (des2.des2 <= tips.num){
        RRF.mat[des2.row.id, 11] <- r4
        t6 <- sqrt(l3*l4)*sqrt(sqrt(l1*l2)+l5)/sqrt(sqrt(l3*l4)+l6)
        RRF.mat[des2.row.id, "t7"] <- t6
      }

    }
  }

  # RRF.mat[is.na(RRF.mat[,10]), 10] <- 0
  # RRF.mat[is.infinite(RRF.mat[,10]), 10] <- 0
  # RRF.mat[is.na(RRF.mat[,11]), 11] <- 0
  # RRF.mat[is.infinite(RRF.mat[,11]), 11] <- 0

  #### adjust relative rates by multiplying the ancestral rate ####
  #### adjust relative times by dividing the ancestral rate ####
  for (nd in (tips.num+1):(tips.num+t$Nnode)){ ## from root to shallowest internal nodes

    row.id <- match(nd, RRF.mat[, 1])
    anc.row <- c(match(nd, RRF.mat[,2]), match(nd, RRF.mat[,3]))

    if (is.na(anc.row[1]) == TRUE && is.na(anc.row[2]) == TRUE){ ## ancestor is root
      r.anc <- 1
      RRF.mat[row.id, c(12,13)] <- RRF.mat[row.id, c(10,11)]*r.anc
      RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc
    }else if (is.na(anc.row[1]) == FALSE && is.na(anc.row[2]) == TRUE){
      r.anc <- RRF.mat[anc.row[1], 12]
      RRF.mat[row.id, c(12,13)] <- RRF.mat[row.id, c(10,11)]*r.anc
      RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc

      if(RRF.mat[row.id, 2] <= tips.num){ ## if one offspring is tip, we need grandparent rate
        grandpa.row <- c(match(RRF.mat[anc.row[1], 1], RRF.mat[, 2]), match(RRF.mat[anc.row[1], 1], RRF.mat[, 3]))
        if (is.na(grandpa.row[1]) == TRUE && is.na(grandpa.row[2]) == TRUE){
          r.anc <- 1
          RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc
        }else{
          r.grandpa <- c(RRF.mat[grandpa.row[1], 12], RRF.mat[grandpa.row[2], 13])
          r.anc <-  r.grandpa[!is.na(r.grandpa)]
          RRF.mat[row.id, 12] <- RRF.mat[row.id, 10]*r.anc
          RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc
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
          RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc
        }
      }
    }else if (is.na(anc.row[1]) == TRUE && is.na(anc.row[2]) == FALSE){
      r.anc <- RRF.mat[anc.row[2], 13]
      RRF.mat[row.id, c(12,13)] <- RRF.mat[row.id, c(10,11)]*r.anc
      RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc

      if(RRF.mat[row.id, 2] <= tips.num){ ## if one offspring is tip, we need grandparent rate
        grandpa.row <- c(match(RRF.mat[anc.row[2], 1], RRF.mat[, 2]), match(RRF.mat[anc.row[2], 1], RRF.mat[, 3]))
        if (is.na(grandpa.row[1]) == TRUE && is.na(grandpa.row[2]) == TRUE){
          r.anc <- 1
          RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc
        }else{
          r.grandpa <- c(RRF.mat[grandpa.row[1], 12], RRF.mat[grandpa.row[2], 13])
          r.anc <-  r.grandpa[!is.na(r.grandpa)]
          RRF.mat[row.id, 12] <- RRF.mat[row.id, 10]*r.anc
          RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc
        }
      }

      if(RRF.mat[row.id, 3] <= tips.num){ ## if one offspring is tip, we need grandparent rate
        grandpa.row <- c(match(RRF.mat[anc.row[2], 1], RRF.mat[, 2]), match(RRF.mat[anc.row[2], 1], RRF.mat[, 3]))
        if (is.na(grandpa.row[1]) == TRUE && is.na(grandpa.row[2]) == TRUE){
          r.anc <- 1
          RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc
        }else{
          r.grandpa <- c(RRF.mat[grandpa.row[1], 12], RRF.mat[grandpa.row[2], 13])
          r.anc <-  r.grandpa[!is.na(r.grandpa)]
          RRF.mat[row.id, 13] <- RRF.mat[row.id, 11]*r.anc
          RRF.mat[row.id, "t7.adjust"] <- RRF.mat[row.id, "t7"]/r.anc
        }
      }
    }
  }

  ## set a rate ratio threshold to avoid extreme rates and node times
  rate.ratio = 20
  node.exceed.ratio.r5 = which(RRF.mat[,"r5.adjust"] > rate.ratio | RRF.mat[,"r5.adjust"] <1/rate.ratio)
  node.exceed.ratio.r6 = which(RRF.mat[,"r6.adjust"] > rate.ratio | RRF.mat[,"r6.adjust"] <1/rate.ratio)
  node.exceed.ratio = union(RRF.mat[node.exceed.ratio.r5, "Des1"], RRF.mat[node.exceed.ratio.r6, "Des2"])

  for (ex in node.exceed.ratio.r5){
    if (RRF.mat[ex, "Des1"] > tips.num){ # not tips
      no.exceed.anc = phangorn::Ancestors(t, RRF.mat[ex, "Des1"], type = "all")
      no.exceed.anc = max(setdiff(no.exceed.anc, node.exceed.ratio))
      RRF.mat[match(RRF.mat[ex, "Des1"], RRF.mat[, "NodeId"]), "t7.adjust"] = RRF.mat[match(no.exceed.anc, RRF.mat[, "NodeId"]), "t7.adjust"]
    }
  }

  for (ex in node.exceed.ratio.r6){
    if (RRF.mat[ex, "Des2"] > tips.num){ # not tips
      no.exceed.anc = phangorn::Ancestors(t, RRF.mat[ex, "Des2"], type = "all")
      no.exceed.anc = max(setdiff(no.exceed.anc, node.exceed.ratio))
      RRF.mat[match(RRF.mat[ex, "Des2"], RRF.mat[, "NodeId"]), "t7.adjust"] = RRF.mat[match(no.exceed.anc, RRF.mat[, "NodeId"]), "t7.adjust"]
    }
  }

  ############ END OF RRF CALCULATION ##############

  ############ BEGIN OF DDBD CALCULATION ##############
  rel.time <- RRF.mat[, "t7.adjust"]/max(RRF.mat[, "t7.adjust"])
  rel.time[rel.time<0] <- 0
  if (anchor.node == 0){
    sf <- 1
  }else{
    sf <- anchor.time/rel.time[match(anchor.node, RRF.mat[, "NodeId"])]
  }

  b.rate.try <- seq(1,10,1)+0.1
  d.rate.try <- seq(1,10,1)
  s.fr.try <-c(0.001,0.01,0.1,0.5,0.9)
  paras.try <- expand.grid(b.rate.try, d.rate.try, s.fr.try)
  paras.try <- paras.try[-which(paras.try$Var1<paras.try$Var2), ] # force birth rate >= death rate

  rel.time.den <- density(rel.time)
  rel.time.den.x <- rel.time.den$x
  rel.time.den.y <- rel.time.den$y
  rel.time.den.x.2 <- rel.time.den.x[rel.time.den.x>=0 & rel.time.den.x<=1]
  rel.time.den.y.2 <- rel.time.den.y[rel.time.den.x>=0 & rel.time.den.x<=1]

  err <- numeric()
  kl.dist <- numeric()

  BD.density <- function(tt, birth.rate, death.rate, rho, root.age=1){
    A <- exp((death.rate - birth.rate)*tt)
    A1 <- exp((death.rate - birth.rate)*root.age)
    Prob.t <- (rho*(birth.rate - death.rate))/(rho*birth.rate + (birth.rate*(1-rho) - death.rate) * A)
    Prob.t1 <- (rho*(birth.rate - death.rate))/(rho*birth.rate + (birth.rate*(1-rho) - death.rate) * A1)
    vt1 <- 1 - 1/rho * Prob.t1*A1
    p1 <- 1/rho * Prob.t^2 * A
    gt <- birth.rate*p1 / vt1

    return(gt)
  }

  for (i in 1:nrow(paras.try)){
    bd.density <- BD.density(tt = rel.time.den.x.2, birth.rate = paras.try[i,1], death.rate = paras.try[i,2], rho = paras.try[i,3], root.age = 1)

    err <- c(err, sqrt(sum((rel.time.den.y.2-bd.density)^2)))
    kl.dist <- c(kl.dist, mean(FNN::KL.divergence(rel.time.den.y.2, bd.density, k=5)))
  }

  err.sort <- sort(err)
  kl.sort <- sort(kl.dist)
  attempt <- 0

  LL.BD <- function(lambda, mu, rho){
    t1 <- 1
    A <- exp((mu-lambda)*rel.time)
    A1 <- exp((mu-lambda)*t1)
    Prob.t <- (rho*(lambda-mu))/(rho*lambda + (lambda*(1-rho) - mu) * A)
    Prob.t1 <- (rho*(lambda-mu))/(rho*lambda + (lambda*(1-rho) - mu) * A1)
    vt1 <- 1 - (1/rho)*Prob.t1*A1
    p1 <- (1/rho)*Prob.t^2*A
    gt <- (lambda*p1)/vt1
    -sum(log(gt))
  }


  if (sampling.frac == 0){
    repeat{
      attempt <- attempt + 1

      if (measure == "KL"){
        paras.start <- paras.try[match(kl.sort[attempt], kl.dist), ]
      }else{
        paras.start <- paras.try[match(err.sort[attempt], err), ]
      }

      names(paras.start) <- c("lambda", "mu", "rho")

      inf.paras <- try(stats4::mle((LL.BD), start = list(lambda = as.numeric(paras.start[1]), mu = as.numeric(paras.start[2]), rho = as.numeric(paras.start[3])),
                                   method = "L-BFGS-B", lower = c(0, 0, 0), upper = c(Inf, Inf, 1)), silent = TRUE)
      if (!(inherits(inf.paras,"try-error")))
        break

      if (attempt >= 50)
        break
    }

    if (is.null(inf.paras) == TRUE){
      stop("Sorry, the best parameter setting cannot be found after trying 50 times.")
    }else{
      b <- inf.paras@coef[1]/sf
      d <- inf.paras@coef[2]/sf
      f <- inf.paras@coef[3]
    }

  }else{
    repeat{
      attempt <- attempt + 1

      if (measure == "KL"){
        paras.start <- paras.try[match(kl.sort[attempt], kl.dist), ]
      }else{
        paras.start <- paras.try[match(err.sort[attempt], err), ]
      }

      names(paras.start) <- c("lambda", "mu", "rho")

      inf.paras <- try(stats4::mle((LL.BD), start = list(lambda = as.numeric(paras.start[1]), mu = as.numeric(paras.start[2]), rho = as.numeric(paras.start[3])),
                                   method = "L-BFGS-B", lower = c(0, 0, 0), upper = c(Inf, Inf, 1)), silent = TRUE)
      if (!(inherits(inf.paras,"try-error")))
        break

      if (attempt >= 50)
        break
    }

    if (is.null(inf.paras) == TRUE){
      stop("Sorry, the best parameter setting cannot be found after trying 50 times.")
    }else{
      b <- inf.paras@coef[1]/sf
      d <- inf.paras@coef[2]/sf
      f <- sampling.frac
    }
  }

  write(paste("birth.rate", "death.rate", "sampling.frac", sep = "\t"), file = paste0(filename, "_ddbd.txt"))
  write(paste(b, d, f, sep = "\t"), file = paste0(filename, "_ddbd.txt"), append = TRUE)

  if (plot == TRUE){ ## node time distribution & curves
    abs.time <- rel.time * sf
    h <- hist(abs.time, breaks = seq(0, 1, 0.025)*sf, col="gray80", freq = FALSE, border = "gray80",
              axes=FALSE, xlab="Node times", ylab="Density", main="", cex.axis = 1.2, cex.lab = 1.2, xaxs="i",yaxs="i")
    axis(1)
    axis(2)

    nt <- seq(0, 1, 0.025)
    bd.density.inf <- BD.density(tt = nt, birth.rate = b*sf, death.rate = d*sf, rho = f, root.age = 1)
    lines(nt*sf, bd.density.inf/sf, col="darkred", lwd = 2)
  }

  return(list("birth.rate" = b, "death.rate" = d, "sampling.fraction" = f))
}

