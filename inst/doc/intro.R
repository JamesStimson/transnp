## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  quiet = TRUE,
  progress = FALSE,
  comment = "#>",
  fig.width=6,
  fig.height=4
)

## ----global_opts, echo=FALSE---------------------------------------------
options(readr.show_progress = FALSE)

## ---- message=FALSE,results="hide"---------------------------------------
library(transnp)
library(ape)
library(phangorn)
library(TransPhylo)
library(treedater)
library(phytools)
library(gridExtra)
library(purrr)
library(ggplot2)
library(visNetwork)

## ------------------------------------------------------------------------
aligns <- list()
aligns <- addAlignment(aligns, "demo_cluster_CL001.fas")
aligns <- addAlignment(aligns, "demo_cluster_CL002.fas")
aligns <- addAlignment(aligns, "demo_cluster_CL005.fas")
aligns <- addAlignment(aligns, "demo_cluster_CL007.fas")

## ------------------------------------------------------------------------
aligns[[1]]

## ---- results='hide'-----------------------------------------------------
basefreqs <- c(0.1719, 0.3286, 0.3276, 0.1719)
set.seed(12345)
mltrees <- createMLTrees(aligns, bf=basefreqs)
# Alternatively, mltrees=lapply(aligns, function(x) estimate.tree(x, bf = basefreqs, maxit=10000,optQ=TRUE, optNni=TRUE))
par(mfrow=c(2,2)); lapply(mltrees,function(tree) {plot(tree); add.scale.bar()}); par(mfrow=c(1,1))
# Improvement to do: coerce plots to increase height and width to show labels better

## ------------------------------------------------------------------------
mltrees[[1]]

## ---- warning = FALSE, messages=FALSE------------------------------------
allDates <- read.csv(system.file("extdata", "demo_dates.csv", package = "transnp", mustWork = TRUE))
sampledates <- setUpDates(aligns, allDates)
timedtrees <- createTimedTrees(mltrees, sampledates, aligns, strictClock = T, meanRateLimits = c(0.3, 0.7))

## ----quiet=TRUE,results="hide"-------------------------------------------
nIter <- 1000
tpTimedTrees <- list()

maxDate <- max(unlist(lapply(sampledates, function(dateList) {return(max(dateList))})))
tpTimedTrees <- lapply(timedtrees, function(tree) {return(ptreeFromPhylo(tree, maxDate))}) 

#listRecords <- infer_multiTTree_shareParam(tpTimedTrees, share=c("neg","off.r","off.p","pi"), mcmcIterations=nIter, startNeg=1, startOff.p = 0.8, startOff.r = 1, startPi=0.9, updateNeg = T, updatePi = F, updateOff.p=FALSE)
listRecords <- simulTransTrees(tpTimedTrees, mean_gen=2, stddev_gen=sqrt(2), mean_sample=2, stddev_sample=sqrt(2), share=c("neg","off.r","off.p","pi"), mcmcIterations=nIter, startNeg=1, startOff.p = 0.8, startOff.r = 1, startPi=0.9, updateNeg = T, updatePi = F, updateOff.p=FALSE)

## ------------------------------------------------------------------------
maptree <- selectTTree(listRecords[[4]],burnin = 0.3)
plotCTree(listRecords[[4]][[maptree]]$ctree)

## ------------------------------------------------------------------------
networkTPlot(listRecords[[4]], mcmcIndex=nIter-10, missLabel="Missed", colours=c('lightblue', 'orange'))

## ------------------------------------------------------------------------
networkTPlot(listRecords[[2]], mcmcIndex=nIter, missLabel="Missed", colours=c('blue', 'orange'))

## ------------------------------------------------------------------------
names(listRecords) = c("Cluster1", "Cluster2", "Cluster5", "Cluster7")
plot_gen_times(listRecords)

## ------------------------------------------------------------------------
plot_times_to_samp(listRecords)

## ------------------------------------------------------------------------
plot_unsampled_cases(listRecords)

## ------------------------------------------------------------------------
gt11=getInfectionTimes(listRecords[[1]][300:1000],7)
name1=listRecords[[1]][[1000]]$ctree$nam[7]
hist(gt11,breaks=15,xlab = paste("estimated time of infection for", name1),main = "")

## ------------------------------------------------------------------------
maptree1=selectTTree(listRecords[[1]],burnin = 0.5)
plotCTree(listRecords[[1]][[maptree1]]$ctree)

## ------------------------------------------------------------------------
plotGenerationTimeDensity(listRecords[[1]], index=2)

## ------------------------------------------------------------------------
plotSampleTimeDensity(listRecords[[1]], index=2)

## ------------------------------------------------------------------------
plotInfectionDateDensity(listRecords[[1]], index=1, overlayDate=2015)

## ------------------------------------------------------------------------
getLikelyInfectors(listRecords[[1]], 'c405')

## ------------------------------------------------------------------------
plotTransTreeSummary(listRecords[[1]], cutoff=0.1, burnin=0.5) 

## ------------------------------------------------------------------------
subRecord1 <- filterTransPair(listRecords[[1]], from=3, to=4)

## ------------------------------------------------------------------------
subRecord2 <- filterInfectedDate(listRecords[[1]], case=1, infDate=2015, keepAfter=TRUE)

## ------------------------------------------------------------------------
testMatrix <- matrix(TRUE, 8, 8)
testMatrix[1,2] <- FALSE
subRecord3 <- filterPossibleWiw(listRecords[[1]], testMatrix)

