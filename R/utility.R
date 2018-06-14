
#' Get generation times
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree
#' @return Vector of times between becoming infected and infecting others (generation times) in the posterior
#' @examples
#' getGenerationTimes(record[[1]]$ctree)
getGenerationTimes <-  function(ctree) { tt=extractTTree(ctree)$ttree;
# 3rd column of ttree lists the infectors; exclude source
infectors=tt[,3]; infectors=infectors[infectors!=0];
# times at which each infector infected others:
infothers=tt[tt[,3]!=0,1];
# times at which each infector was herself infected:
gotinfd=tt[infectors,1];
return(infothers-gotinfd);}

#' Get the number of unsampled cases
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree
#' @return The number of unsampled cases in ctree
#' @examples
#' getNumberUnsampled(ctree)
getNumberUnsampled <- function(ctree) {tt=extractTTree(ctree); return(sum(is.na(tt$ttree[,2])))}

#' Create list of wiw information in order to compute transmission tree distances
#' @param record  MCMC output produced by inferTTree
#' @param skipnum Number of record entries to skip, ie skipnum=10 uses the 1st, 11th, 21st ... elements of the record
#' @return list of MRCI information required by wiwTreeDist, one entry for each transmission tree that is included
getTTreeDistInfo <- function(record,skipnum=1) {
  ind=seq(from=1, by=skipnum,to=length(record))
  record=record[ind]
  matList <- lapply(1:length(record), function(x) {
    info <- extractTTree(record[[x]]$ctree)$ttree
    wiw <- cbind(info[,3],1:length(info[,1]))
    findMRCIs(wiw)$mrciDepths
  })
  return(matList)
}

#' Get infection likelihoods for a given case
#' @param record  MCMC output produced by inferTTree
#' @param case  Index or label can be supplied
#' @param burnin Proportion of record entries to use
#' @export
getLikelyInfectors <- function(record, case, burnin=0.5){
  # Who-infected-who matrix from TransPhylo
  myMat <- computeMatWIW(record, burnin)
  return(myMat[,case])
}

#' Get infection-to-sampling times
#' @param ctree TransPhylo tree
#' @return A vector of posterior times to sampling
#' @examples
#' getTimesToSampling(ctree)
getTimesToSampling <-  function(ctree) { tt=extractTTree(ctree)$ttree;
ns=sum(!is.na(tt[,2]))
return(tt[1:ns,2]-tt[1:ns,1])
}

#' Get infection dates
#' @param ctree TransPhylo tree
#' @return A vector of posterior infection dates
#' @examples
#' getInfectionDates(ctree)
getInfectionDates <-  function(ctree) { tt=extractTTree(ctree)$ttree;
ns=sum(!is.na(tt[,2]))
return(tt[1:ns,1])
}

#' add alignments to list
#' @param alist (possibly empty) list of alignments in DNAbin format
#' @param fileName fasta file name
#' @return list of alignments
#' @export
addAlignment <- function(alist, fileName){
  alignment <- read.dna(system.file("extdata", fileName, package = "transnp", mustWork = TRUE), format="fasta")
  alist[[length(alist)+1]] <- alignment
  return (alist)
}

#' create list of dates to match alignments
#' @param alist (possibly empty) list of alignments in DNAbin format
#' @param allDates data frame containg all IDs and dates as read from a csv file
#' @return list of dates
#' @export
setUpDates <- function(alist, allDates){
  dateList <- list()
  for (align in alist){
    theseDates <- numeric(length(labels(align)))
    for (i in seq(length(labels(align)))){
      theseDates[[i]] <- allDates[match(labels(align)[[i]], allDates[,1]), 2]
    }
    dateList[[length(dateList)+1]] <- theseDates
  }
  return(dateList)
}





