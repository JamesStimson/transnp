#' Get infection times
#' @param record MCMC output produced by inferTTree
#' @param k Case whose posterior infection times are to be extracted. Either an integer or a string matching one of the case names in the data
#' @return A vector of posterior infection times for case k
#' @export
#' @examples
#' getInfectionTimes(record)
getInfectionTimes <- function(record,k) {
  if (is.numeric(k)) {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); return(tt$ttree[k,1])},FUN.VALUE=1);
    return(mytimes)}
  else {
    mytimes= vapply(1:length(record), function(x) { tt=extractTTree(record[[x]]$ctree); ii=which(tt$nam==k); return(tt$ttree[ii,1])},FUN.VALUE=1);
    return(mytimes)}
}


#' Filter record on specified directed transmission pair
#' @param record MCMC output produced by inferTTree
#' @param from From case index, either an integer or a string matching one of the case names in the data
#' @param to To case index, either an integer or a string matching one of the case names in the data
#' @return Filtered record
#' @export
filterTransPair <- function(record, from, to){
  rList <- list()
  count <- 0
  for  (entry in record){
    tt <- extractTTree(entry$ctree)
    if(!is.numeric(from)) from <- which(tt$nam==from)
    if(!is.numeric(to)) to <- which(tt$nam==to)
    if (tt$ttree[to, 3] == from){
      count <- count + 1
      rList[[count]] <- entry
    }
  }
  return(rList)
}

#' Filter record on specified date of infection
#' @param record MCMC output produced by inferTTree
#' @param case Case index, either an integer or a string matching one of the case names in the data
#' @param infDate Date for filtering
#' @param keepAfter Flag, TRUE to keep cases infected after infDate, FALSE to keep cases infected on or before
#' @return Filtered record
#' @export
filterInfectedDate <- function(record, case, infDate, keepAfter=TRUE){
  rList <- list()
  count <- 0
  for  (entry in record){
    tt <- extractTTree(entry$ctree)
    if(!is.numeric(case)) case <- which(tt$nam==case)
    if ((keepAfter && (tt$ttree[case, 1] > infDate)) || (!keepAfter && (tt$ttree[case, 1] <= infDate))){
      count <- count + 1
      rList[[count]] <- entry
    }
  }
  return(rList)
}

#' Filter record on specified allowable transmission matrix
#' @param record MCMC output produced by inferTTree
#' @param pMatrix Matrix of allowable transmissions, TRUE for allowed, FALSE for not allowed
#' @return Filtered record
#' @export
filterPossibleWiw <- function(record, pMatrix){
  rList <- list()
  count <- 0
  for  (entry in record){
    tt <- extractTTree(entry$ctree)
    allow <- TRUE
    for (i in seq(1:length(tt$nam))){
      infector <- tt$ttree[i,3]
      if (infector > length(tt$nam)) next
      if (!pMatrix[tt$ttree[i,3], i]){
        allow <- FALSE
        break
      }
    }
    if (allow){
      count <- count + 1
      rList[[count]] <- entry
    }
  }
  return(rList)
}
