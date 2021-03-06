# Note latest github TransPhylo needed
# install_github("xavierdidelot/TransPhylo")

#' create maximum likelihood trees from a list of alignments
#' @param alist list of alignments in DNAbin format
#' @param maxit Max number of iterations in optimisation
#' @param optNni (TRUE) optimise the topology of the tree
#' @param ... further parameters for optim.pml
#' @param untangle flag as to whether to untangle tree, may be needed for small trees
#' @return list of maximum likelihood trees in phylo format
#' @export
createMLTrees <- function(alist, maxit=10000, optNni=TRUE, bf=basefreqs, untangle=TRUE, ...) {
  returnTrees <- list()
  #set.seed(12345)
  for (a in alist){
    mlTree <- estimatetree(a, bf=basefreqs)
    # adjust to make branch lengths ~ number SNPs
    mlTree$edge.length <- mlTree$edge.length * ncol(a)
    # untangle problematic trees
    if(untangle){
      mlTree <- untangle(mlTree)
    }
    # fix any very small, negative or zero edge lengths
    mlTree$edge.length[which(mlTree$edge.length < 1e-5)]<- 1e-5
    returnTrees[[length(returnTrees)+1]] <- mlTree
  }
  return (returnTrees)
}

#' create timed trees from a list of maximum likelihood trees
#' @param mltrees list of trees as created by createMLTrees
#' @param sampledates vector of sample dates
#' @param alist list of alignments in DNAbin format
#' @param (TRUE) set whether to use strict clock in treedater
#' @param meanRateLimits vector of length 2 specifying upper and lower rate limits for treedater
#' @return list of timed trees in phylo format
#' @export
createTimedTrees <- function(mltrees, sampledates, alist, strictClock = T, meanRateLimits = c(0.3, 0.7)) {
  timedTrees <- list()
  for (i in seq(length(mltrees))){
    theLabels <- labels(alist[[i]])
    theDates <- sampledates[[i]]
    names(theDates) <- theLabels
    newTree <- dater(tre=mltrees[[i]], sts=theDates, s=ncol(alist[[i]]), strictClock = T, meanRateLimits = c(0.3, 0.7))
    # fix any very small, negative or zero edge lengths
    newTree$edge.length[which(newTree$edge.length < 1e-5)]<- 1e-5

    timedTrees[[length(timedTrees)+1]] <- newTree
  }
  return (timedTrees)
}

#' estimate a tree
#' @param aln Alignment in DNAbin format
#' @param maxit Max number of iterations in optimisation
#' @param optNni (TRUE) optimise the topology of the tree
#' @param ... further parameters for optim.pml
#' @return tree object of type phylo
#' @examples test=estimatetree(aln)
#' myModel <- setDatesFromFile(myModel, datefile)
estimatetree <- function(aln, maxit=10000, optNni=TRUE, bf=basefreqs, ...) {
  # create a neighbour-joining tree with the default K80 model to start
  dd <- dist.dna(aln, model="N")
  tre.ini <- bionj(dd)
  # fix any negative or zero edges
  tre.ini$edge.length[which(tre.ini$edge.length <= 0)]<- 1e-5
  # initial likelihood
  fit.ini <- pml(tre.ini,data=as.phyDat(aln), bf=bf, model="K80", k=4)
  # compute ML tree
  fit <- optim.pml(fit.ini, maxit=maxit, optNni=optNni, ...)
  return(fit$tree)
}

#' estimate a tree
#' @param aln Alignment in DNAbin format
#' @param maxit Max number of iterations in optimisation
#' @param optNni (TRUE) optimise the topology of the tree
#' @param ... further parameters for optim.pml
#'
estimate.tree <- function(aln,maxit=10000,addCols=TRUE, optNni=TRUE, bf=basefreqs, ...) {
  # create a neighbour-joining tree with the default K80 model to start
  # pml does not work with tiny alignments. the apparent distances are too big (too high a portion of the tiny "genome")
  # this is a bit silly, but it seems to work:
  if (addCols==TRUE) {
    fakerow=sample(c("a","c","g","t"), size=5000, replace=TRUE, prob=basefreqs)
    fakedata=matrix(rep(fakerow,nrow(aln)),nrow = nrow(aln),ncol=length(fakerow),byrow = TRUE)
    aln=as.character(aln)
    aln=as.DNAbin(cbind(aln, fakedata))
  }
  dd=dist.dna(aln,model="raw")
  tre.ini=bionj(dd)
  # fix any negative or zero edges
  tre.ini$edge.length[which(tre.ini$edge.length <= 0)]<- 1e-5
  # initial likelihood
  fit.ini <- pml(tre.ini,data=as.phyDat(aln),bf= bf, model="K80", k=4)
  # compute ML tree
  fit <- optim.pml(fit.ini, maxit=maxit, optNni=optNni, ...)
  tree=fit$tree
  tree$edge.length=tree$edge.length*ncol(aln)
  return(tree)
}

#' Simultaneously infer transmission trees given phylogenetic trees constructed from clusters of sequences.
#' Wrapper for infer_multiTTree_shareParam, taking mean and std deviations rather than shape and scale for the Gamma parameters
#' @param ptree_lst List of phylogenetic tree
#' @param mean_gen Mean of the Gamma probability density function representing the generation time
#' @param stddev_gen Standard deviation of the Gamma probability density function representing the generation time
#' @param mean_sample Mean of the Gamma probability density function representing the sampling time
#' @param stddev_sample Standard deviation of the Gamma probability density function representing the sampling time
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param share Character vector of parameters to be shared. For example, share = c("off.r", "off.p") would
#' share the offspring distribution. Allowed parameter names are "neg", "off.r", "off.p" and "pi".
#' @param startCTree_lst Optional combined list of trees to start from
#' @param updateTTree Whether or not to update the transmission tree
#' @param optiStart Whether or not to optimise the MCMC start point
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return posterior sample set of transmission trees for all clusters
#' @export
simulTransTrees = function(timedTrees, mean_gen=1, stddev_gen=1, mean_sample=1, stddev_sample=1,
                           share=c("neg","off.r","off.p","pi"), mcmcIterations=nIter, startNeg=1, startOff.p = 0.8, startOff.r = 1, startPi=0.9, updateNeg = T, updatePi = F, updateOff.p=FALSE){

  # mean is shape*scale, st dev is sqrt(shape)*scale, so:
  shape_gen <- (mean_gen/stddev_gen)^2
  scale_gen <- (stddev_gen^2)/mean_gen
  shape_sample <- (mean_sample/stddev_sample)^2
  scale_sample <- (stddev_sample^2)/mean_sample

  record <- infer_multiTTree_shareParam(timedTrees, share=share, w.shape=shape_gen, w.scale=scale_gen, ws.shape=shape_sample, ws.scale=scale_sample,
                                        mcmcIterations=mcmcIterations, startNeg=startNeg, startOff.p = startOff.p, startOff.r = startOff.r, startPi=startPi, updateNeg=updateNeg, updatePi=updatePi, updateOff.p=updateOff.p)
  return(record)
}

#' Simultaneously infer transmission trees given phylogenetic trees constructed from clusters of sequences.
#' User can specify any parameter(s) to be shared by providing a character vector of parameter names to
#' the argument "share".
#' @param ptree_lst List of phylogenetic tree
#' @param w.shape Shape parameter of the Gamma probability density function representing the generation time
#' @param w.scale Scale parameter of the Gamma probability density function representing the generation time
#' @param ws.shape Shape parameter of the Gamma probability density function representing the sampling time
#' @param ws.scale Scale parameter of the Gamma probability density function representing the sampling time
#' @param mcmcIterations Number of MCMC iterations to run the algorithm for
#' @param thinning MCMC thinning interval between two sampled iterations
#' @param startNeg Starting value of within-host coalescent parameter Ne*g
#' @param startOff.r Starting value of parameter off.r
#' @param startOff.p Starting value of parameter off.p
#' @param startPi Starting value of sampling proportion pi
#' @param updateNeg Whether of not to update the parameter Ne*g
#' @param updateOff.r Whether or not to update the parameter off.r
#' @param updateOff.p Whether or not to update the parameter off.p
#' @param updatePi Whether or not to update the parameter pi
#' @param share Character vector of parameters to be shared. For example, share = c("off.r", "off.p") would
#' share the offspring distribution. Allowed parameter names are "neg", "off.r", "off.p" and "pi".
#' @param startCTree_lst Optional combined list of trees to start from
#' @param updateTTree Whether or not to update the transmission tree
#' @param optiStart Whether or not to optimise the MCMC start point
#' @param dateT Date when process stops (this can be Inf for fully simulated outbreaks)
#' @return posterior sample set of transmission trees for all clusters
#' @author Yuanwei Xu
infer_multiTTree_shareParam = function(ptree_lst,w.shape=2,w.scale=1,ws.shape=w.shape,ws.scale=w.scale,mcmcIterations=1000,
                                       thinning=1,startNeg=100/365,startOff.r=1,startOff.p=0.5,startPi=0.5,
                                       updateNeg=TRUE,updateOff.r=TRUE,updateOff.p=FALSE,updatePi=TRUE,
                                       share=NULL,
                                       startCTree_lst=rep(NA,length(ptree_lst)),updateTTree=TRUE,optiStart=TRUE,dateT=Inf) {

  ptree_lst <- purrr::map(ptree_lst, function(x) within(x, ptree[,1] <- ptree[,1]+runif(nrow(ptree))*1e-10))
  #MCMC algorithm
  neg <- startNeg
  off.r <- startOff.r
  off.p <- startOff.p
  pi <- startPi

  ctree_lst <- vector("list", length(ptree_lst))
  for(k in seq_along(ptree_lst)){ # starting ctree
    if(is.na(startCTree_lst[[k]]))
      ctree_lst[[k]] <- makeCtreeFromPTree(ptree_lst[[k]],ifelse(optiStart,off.r,NA),off.p,neg,pi,w.shape,w.scale,ws.shape,ws.scale,dateT)
    else
      ctree_lst[[k]] <- startCTree_lst[[k]]
  }
  ntree <- length(ptree_lst)
  neg_lst <- as.list(rep(neg, ntree))
  off.r_lst <- as.list(rep(off.r, ntree))
  off.p_lst <- as.list(rep(off.p, ntree))
  pi_lst <- as.list(rep(pi, ntree))
  not_share <- setdiff(c("neg", "off.r", "off.p", "pi"), share)

  one_update <- function(ctree, pTTree, pPTree, neg, off.r, off.p, pi, not_share){
    # Get a copy of current ttree
    ttree <- extractTTree(ctree)

    if (updateTTree) {
      #Metropolis update for transmission tree
      prop <- .proposal(ctree$ctree)
      ctree2 <- list(ctree=prop$tree,nam=ctree$nam)
      ttree2 <- extractTTree(ctree2)
      pTTree2 <- probTTree(ttree2$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT)
      pPTree2 <- probPTreeGivenTTree(ctree2,neg)
      if (log(runif(1)) < log(prop$qr)+pTTree2 + pPTree2-pTTree-pPTree)  {
        ctree <- ctree2
        ttree <- ttree2
        pTTree <- pTTree2
        pPTree <- pPTree2
      }
    }

    if (("neg" %in% not_share) && updateNeg) {
      #Metropolis update for Ne*g, assuming Exp(1) prior
      neg2 <- abs(neg + (runif(1)-0.5)*0.5)
      pPTree2 <- probPTreeGivenTTree(ctree,neg2)
      if (log(runif(1)) < pPTree2-pPTree-neg2+neg)  {neg <- neg2;pPTree <- pPTree2}
    }

    if (("off.r" %in% not_share) && updateOff.r) {
      #Metropolis update for off.r, assuming Exp(1) prior
      off.r2 <- abs(off.r + (runif(1)-0.5)*0.5)
      pTTree2 <- probTTree(ttree$ttree,off.r2,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT)
      if (log(runif(1)) < pTTree2-pTTree-off.r2+off.r)  {off.r <- off.r2;pTTree <- pTTree2}
    }

    if (("off.p" %in% not_share) && updateOff.p) {
      #Metropolis update for off.p, assuming Unif(0,1) prior
      off.p2 <- abs(off.p + (runif(1)-0.5)*0.1)
      if (off.p2>1) off.p2=2-off.p2
      pTTree2 <- probTTree(ttree$ttree,off.r,off.p2,pi,w.shape,w.scale,ws.shape,ws.scale,dateT)
      if (log(runif(1)) < pTTree2-pTTree)  {off.p <- off.p2;pTTree <- pTTree2}
    }

    if (("pi" %in% not_share) && updatePi) {
      #Metropolis update for pi, assuming Unif(0.01,1) prior
      pi2 <- pi + (runif(1)-0.5)*0.1
      if (pi2<0.01) pi2=0.02-pi2
      if (pi2>1) pi2=2-pi2
      pTTree2 <- probTTree(ttree$ttree,off.r,off.p,pi2,w.shape,w.scale,ws.shape,ws.scale,dateT)
      if (log(runif(1)) < pTTree2-pTTree)  {pi <- pi2;pTTree <- pTTree2}
    }

    list(ctree=ctree, pTTree=pTTree, pPTree=pPTree, neg=neg, off.r=off.r, off.p=off.p, pi=pi)
  }

  one_update_share <- function(ctree_lst, pTTree_lst, pPTree_lst, neg_lst, off.r_lst, off.p_lst, pi_lst, share){
    ttree_lst <- purrr::map(ctree_lst, extractTTree)
    ttree_lst <- purrr::map(ttree_lst, "ttree") # list of ttree matrices

    if(("neg" %in% share) && updateNeg){
      neg <- neg_lst[[1]]
      #Metropolis update for Ne*g, assuming Exp(1) prior
      neg2 <- abs(neg + (runif(1)-0.5)*0.5)
      pPTree <- purrr::flatten_dbl(pPTree_lst)
      pPTree2 <- purrr::map_dbl(ctree_lst, probPTreeGivenTTree, neg = neg2)
      if (log(runif(1)) < sum(pPTree2)-sum(pPTree)-neg2+neg){
        neg_lst <- as.list(rep(neg2, ntree))
        pPTree_lst <- as.list(pPTree2)
      }
    }

    if(("off.r" %in% share) && updateOff.r) {
      off.r <- off.r_lst[[1]]
      #Metropolis update for off.r, assuming Exp(1) prior
      off.r2 <- abs(off.r + (runif(1)-0.5)*0.5)
      pTTree <- purrr::flatten_dbl(pTTree_lst)
      pTTree2 <- purrr::pmap_dbl(list(ttree=ttree_lst, pOff=off.p_lst, pi=pi_lst), probTTree, rOff=off.r2,
                                 shGen=w.shape, scGen=w.scale, shSam=ws.shape, scSam=ws.scale, dateT=dateT)
      if (log(runif(1)) < sum(pTTree2)-sum(pTTree)-off.r2+off.r){
        off.r_lst <- as.list(rep(off.r2, ntree))
        pTTree_lst <- as.list(pTTree2)
      }
    }

    if(("off.p" %in% share) && updateOff.p) {
      off.p <- off.p_lst[[1]]
      #Metropolis update for off.p, assuming Unif(0,1) prior
      off.p2 <- abs(off.p + (runif(1)-0.5)*0.1)
      if (off.p2>1) off.p2=2-off.p2
      pTTree <- purrr::flatten_dbl(pTTree_lst)
      pTTree2 <- purrr::pmap_dbl(list(ttree=ttree_lst, rOff=off.r_lst, pi=pi_lst), probTTree, pOff=off.p2,
                                 shGen=w.shape, scGen=w.scale, shSam=ws.shape, scSam=ws.scale, dateT=dateT)
      if (log(runif(1)) < sum(pTTree2)-sum(pTTree)){
        off.p_lst <- as.list(rep(off.p2, ntree))
        pTTree_lst <- as.list(pTTree2)
      }
    }

    if(("pi" %in% share) && updatePi){
      pi <- pi_lst[[1]]
      #Metropolis update for pi, assuming Unif(0.01,1) prior
      pi2 <- pi + (runif(1)-0.5)*0.1
      if (pi2<0.01) pi2=0.02-pi2
      if (pi2>1) pi2=2-pi2
      pTTree <- purrr::flatten_dbl(pTTree_lst)
      pTTree2 <- purrr::pmap_dbl(list(ttree=ttree_lst, rOff=off.r_lst, pOff=off.p_lst), probTTree, pi=pi2,
                                 shGen=w.shape, scGen=w.scale, shSam=ws.shape, scSam=ws.scale, dateT=dateT)
      if (log(runif(1)) < sum(pTTree2)-sum(pTTree)){
        pi_lst <- as.list(rep(pi2, ntree))
        pTTree_lst <- as.list(pTTree2)
      }
    }

    list(ctree=ctree_lst, pTTree=pTTree_lst, pPTree=pPTree_lst, neg=neg_lst, off.r=off.r_lst, off.p=off.p_lst, pi=pi_lst)
  }

  # Create outpout data structure: nested list
  record <- vector('list', ntree)
  for(i in seq_along(record)){
    record[[i]] <- vector("list", mcmcIterations/thinning)
  }

  # Initialize MCMC state in multi-tree space
  pTTree_lst <- vector("list", ntree)
  pPTree_lst <- vector("list", ntree)
  for(k in 1:ntree){
    ttree <- extractTTree(ctree_lst[[k]])
    pTTree_lst[[k]] <- probTTree(ttree$ttree,off.r,off.p,pi,w.shape,w.scale,ws.shape,ws.scale,dateT)
    pPTree_lst[[k]] <- probPTreeGivenTTree(ctree_lst[[k]],neg)
  }
  mcmc_state <- list(ctree=ctree_lst, pTTree=pTTree_lst, pPTree=pPTree_lst,
                     neg=neg_lst, off.r=off.r_lst, off.p=off.p_lst, pi=pi_lst)

  #Main MCMC loop
  pb <- txtProgressBar(min=0,max=mcmcIterations,style = 3)
  for (i in 1:mcmcIterations) {

    # Update shared parameters
    out <- with(mcmc_state, one_update_share(ctree, pTTree, pPTree, neg, off.r, off.p, pi, share))
    mcmc_state[["ctree"]] <- out[["ctree"]]
    mcmc_state[["pTTree"]] <- out[["pTTree"]]
    mcmc_state[["pPTree"]] <- out[["pPTree"]]
    mcmc_state[["neg"]] <- out[["neg"]]
    mcmc_state[["off.r"]] <- out[["off.r"]]
    mcmc_state[["off.p"]] <- out[["off.p"]]
    mcmc_state[["pi"]] <- out[["pi"]]

    # Update unshared parameters
    state_new <- purrr::pmap(mcmc_state, one_update, not_share=not_share)

    if (i%%thinning == 0) {
      #Record things
      setTxtProgressBar(pb, i)
      for(k in seq_along(ctree_lst)){
        record[[k]][[i/thinning]]$ctree <- state_new[[k]]$ctree
        record[[k]][[i/thinning]]$pTTree <- state_new[[k]]$pTTree
        record[[k]][[i/thinning]]$pPTree <- state_new[[k]]$pPTree
        record[[k]][[i/thinning]]$neg <- state_new[[k]]$neg
        record[[k]][[i/thinning]]$off.r <- state_new[[k]]$off.r
        record[[k]][[i/thinning]]$off.p <- state_new[[k]]$off.p
        record[[k]][[i/thinning]]$pi <- state_new[[k]]$pi
        record[[k]][[i/thinning]]$w.shape <- w.shape
        record[[k]][[i/thinning]]$w.scale <- w.scale
        record[[k]][[i/thinning]]$ws.shape <- ws.shape
        record[[k]][[i/thinning]]$ws.scale <- ws.scale
        record[[k]][[i/thinning]]$source <- with(state_new[[k]]$ctree, ctree[ctree[which(ctree[,4]==0),2],4])
        if (record[[k]][[i/thinning]]$source<=length(state_new[[k]]$ctree$nam))
          record[[k]][[i/thinning]]$source=state_new[[k]]$ctree$nam[record[[k]][[i/thinning]]$source]
        else record[[k]][[i/thinning]]$source='Unsampled'
      }
    }
    # Assign updated state to current state
    mcmc_state <- purrr::transpose(state_new)

  }#End of main MCMC loop

  #close(pb)
  return(record)
}

#' Build vector 'host' indicating in which host each node of the ctree is found
.computeHost = function(ctree)  {
  fathers <- rep(NA, nrow(ctree))
  fathers[ctree[ ,2] + 1] <- 1:nrow(ctree)
  fathers[ctree[ ,3] + 1] <- 1:nrow(ctree)
  fathers <- fathers[-1]
  host <- rep(0, nrow(ctree))
  nsam <- sum(ctree[ ,2] == 0&ctree[ ,3] == 0)
  for (i in 1:nsam) {
    j <- i
    while (1)  {
      if (host[j]>0) print('Error: two leaves in same host')
      host[j] <- i
      j <- fathers[j]
      if (ctree[j,3] == 0) break
    }
  }
  if (nsam==1) return(host)

  dispo=nsam+1
  f <- which( ctree[,3] == 0 & ctree[,2]>0 & (1:nrow(ctree))<nrow(ctree) ) #transmission events other than root
  for (i in 1:length(f)) {
    j <- f[i]
    tocol <- c()
    while (ctree[fathers[j],3]>0&&fathers[j]<nrow(ctree))  {
      tocol <- c(tocol,j)
      j <- fathers[j]
    }
    if (host[j]==0) {host[j]=dispo;dispo=dispo+1}
    host[tocol] <- host[j]
  }
  return(host)
}
