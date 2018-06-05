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
#' @param from from case index
#' @param to to case index
#' @return Filtered record
filterTransPair <- function(record, from, to){
  rList <- list()
  count = 0
  for  (entry in record){
    tt <- extractTTree(entry$ctree)
    if (tt$ttree[to,3] == from){
      count = count + 1
      rList[[count]] <- entry
    }
  }
  return(rList)
}
