# 2018 version for transnp (was SNPphylo) package


#' Get generation times
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
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
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
#' @return The number of unsampled cases in ctree
#' @examples
#' getNumberUnsampled(ctree)
getNumberUnsampled <- function(ctree) {tt=extractTTree(ctree); return(sum(is.na(tt$ttree[,2])))}

#' Simple visualisation of transmission tree with visNetwork
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree, for example record[[1]]$ctree where record is the output of inferTTree
networkPlot <- function(ctree) {
  tt=extractTTree(ctree)
  info1=tt$ttree
  SimpleWiw <- cbind(info1[,3],1:length(info1[,1])) # infector, infectee
  nodes <- data.frame(id = 1:nrow(info1), label=c(tt$nam,(length(tt$nam)+1) : nrow(info1)))
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to")
  visNetwork(nodes, edges)
}


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


#' dynamic plot of the transmission network with edges weighted according  to likelihood of transmission
#' @param thisRecord Posterior sample set of transmission trees for all clusters
#' @param mcmcIndex Sample set index of network to be displayed
#' @param missLabel Label to be used for missing cases
#' @param colours Colours for sampled and unsampled nodes, vector of 2 colours
#' @return NULL
#' @export
#' @examples
#' networkTPlot(record)
networkTPlot <- function(thisRecord, mcmcIndex=1, missLabel="Unsampled", colours=c('lightblue', 'orange')) {
  #ctree <- thisRecord[[round(length(thisRecord)/2)]]$ctree
  ctree <- thisRecord[[mcmcIndex]]$ctree
  tt <- extractTTree(ctree)
  info1 <- tt$ttree
  myMat <- computeMatWIW(thisRecord, burnin=0.5)

  SimpleWiw <- cbind(info1[,3],1:length(info1[,1])) # infector, infectee
  nodes <- data.frame(id = 1:nrow(info1), label=c(tt$nam, rep(missLabel,nrow(info1)-length(tt$nam))),
                      color=c(rep(colours[[1]],length(tt$nam)),rep(colours[[2]],nrow(info1)-length(tt$nam))),
                      group=c(rep("Sampled",length(tt$nam)),rep(missLabel,nrow(info1)-length(tt$nam))))
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to")
  width <- rep(0.4, nrow(edges))
  for (w in seq(length(width))){
    width[w] <- 10*getMat(edges$from[[w]], edges$to[[w]], myMat)
  }
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to", width=width)

  if(nrow(info1)!=length(tt$nam)){
  visNetwork(nodes, edges) %>%
    visGroups(groupname = "Sampled", color = "lightblue") %>%
    visGroups(groupname = "Unsampled", color = "orange") %>%
    #visLegend(width = 0.1, position = "left", main = "Legend")%>%
    visNodes(font = list(size = 24))
  }
  else{# doesnt work well with just one group!
    visNetwork(nodes, edges)%>%
      visNodes(font = list(size = 24))#%>%
      #visEvents(doubleClick = "itplot_dblclick") # THIS MAKES THE IMAGE DISAPPEAR!
  }
}

#' Plot histogram of generation times for each cluster
#' @param cls_record Posterior sample set of transmission trees for all clusters
#' @param ... parameters to be passed to grid.arrange
#' @return NULL
#' @export
plot_gen_times <- function(cls_record, ...){
  gs <- list()
  for(nm in names(cls_record)){
    gen_times <- map(cls_record[[nm]], ~ getGenerationTimes(.$ctree))
    gs[[nm]] <- ggplot(data.frame(gen_times = unlist(gen_times)), aes(gen_times)) +
      geom_histogram(bins = 20) + ggtitle(nm)

    # The use of data.frame() in ggplot creates a new plot env for each ggplot object
    # See here https://stackoverflow.com/questions/39799886/r-assigning-ggplot-objects-to-list-in-loop/39800861#39800861
    # for explanation of why the following doesn't work properly
    # https://stackoverflow.com/questions/49978925/saving-ggplot-object)
    #gs[[nm]] <- qplot(unlist(gen_times), geom = "histogram", main = nm)
  }
  gridExtra::grid.arrange(grobs = gs, ...)
}

#' Get infection times
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

plotInfectionDateDensity <- function(singleRecord, index, overlayDate){
  dateList <- lapply(1:length(singleRecord), function(x) {
    info <- getInfectionDates(singleRecord[[x]]$ctree)[[index]]
    return(info)
  })
  d <- density(as.numeric(dateList))
  plot(d)
  abline(v=overlayDate, col='red', lwd=2)
}

#' Plot histogram of times to sampling for each cluster
#' @param cls_record Posterior sample set of transmission trees for all clusters
#' @param ... parameters to be passed to grid.arrange
#' @return NULL
#' @export
plot_times_to_samp <- function(cls_record, ...){
  gs <- list()
  for(nm in names(cls_record)){
    times_to_samp <- map(cls_record[[nm]], ~ getTimesToSampling(.$ctree))
    gs[[nm]] <- ggplot(data.frame(times_to_samp = unlist(times_to_samp)), aes(times_to_samp)) +
      geom_histogram(bins = 20) + ggtitle(nm)
  }
  gridExtra::grid.arrange(grobs = gs, ...)
}

#' Plot histogram of the number of unsampled cases for each cluster
#' @param cls_record Posterior sample set of transmission trees for all clusters
#' @return NULL
#' @export
plot_unsampled_cases <- function(cls_record, ...){
  gs <- list()
  for(nm in names(cls_record)){
    unsampled_cases <- map(cls_record[[nm]], ~ getNumberUnsampled(.$ctree))
    gs[[nm]] <- ggplot(data.frame(unsampled_cases = unlist(unsampled_cases)), aes(unsampled_cases)) +
      geom_histogram(bins = 20) + ggtitle(nm)
  }
  gridExtra::grid.arrange(grobs = gs, ...)
}




