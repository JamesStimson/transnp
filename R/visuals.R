
#' Simple visualisation of transmission tree with visNetwork
#' @param ctree Combined tree object 'ctree' from the record produced by inferTTree
networkPlot <- function(ctree) {
  tt=extractTTree(ctree)
  info1=tt$ttree
  SimpleWiw <- cbind(info1[,3],1:length(info1[,1])) # infector, infectee
  nodes <- data.frame(id = 1:nrow(info1), label=c(tt$nam,(length(tt$nam)+1) : nrow(info1)))
  edges <- data.frame(from = SimpleWiw[,1], to = SimpleWiw[,2],
                      arrows="to")
  visNetwork(nodes, edges)
}


#' Plot transmission network with posterior probabilities at or above given cutoff
#' @param record  Posterior MCMC sample set from TransPhylo
#' @param cutoff  Minimum probabilty for edge to be included
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param nodeColour Colour for all nodes
#' @param fontSize Font side for node labels
#' @export
plotTransTreeSummary <- function(record, cutoff=1, burnin=0.5, nodeColour='lightblue', fontSize=24) {
  tt <- extractTTree(record[[1]]$ctree)
  numSampled <- length(tt$nam)
  numElements <- length(record)

  # Who-infected-who matrix from TransPhylo
  myMat <- computeMatWIW(record, burnin)
  # Anything less than the cutoff is set to zero
  myMat[myMat < cutoff] <- 0.0

  from <- c()
  to <- c()
  eLabels <- c()

  nodes <- data.frame(id=1:numSampled, label=tt$nam)
  for(i in 1:nrow(myMat)) {
    for(j in 1:ncol(myMat)) {
      if (myMat[i,j] > 0.0){
        from <- c(from, i)
        to <- c(to, j)
        eLabels <- c(eLabels, round(myMat[i,j], 3))
      }
    }
  }
  width <- rep(1, length(from))
  for (w in seq(length(width))){
    width[w] <- 10*getMat(from[[w]], to[[w]], myMat)
  }
  edges <- data.frame(from = from, to = to, arrows="to", width=width, title=eLabels)
  # Do the plot
  visNetwork(nodes, edges)%>%
    visNodes(color=nodeColour, font=list(size=fontSize))
}

#' Dynamic plot of the transmission network with edges weighted according to likelihood of transmission
#' @param thisRecord Posterior MCMC sample set from TransPhylo
#' @param mcmcIndex Sample set index of network to be displayed
#' @param burnin Proportion of the MCMC output for edge weights to be discarded as burnin
#' @param missLabel Label to be used for missing/unsampled cases
#' @param colours Colours for sampled and unsampled nodes, vector of 2 colours
#' @return NULL
#' @export
#' @examples
#' networkTPlot(record[[5]])
networkTPlot <- function(thisRecord, mcmcIndex=1, burnin=0.5, missLabel="Unsampled", colours=c('lightblue', 'orange'), fontSize=24) {
  #ctree <- thisRecord[[round(length(thisRecord)/2)]]$ctree
  ctree <- thisRecord[[mcmcIndex]]$ctree
  tt <- extractTTree(ctree)
  info1 <- tt$ttree
  myMat <- computeMatWIW(thisRecord, burnin)

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
    visNodes(font = list(size = fontSize))
  }
  else{
    visNetwork(nodes, edges)%>%
      visNodes(font = list(size = fontSize))#%>%
      #visEvents(doubleClick = "itplot_dblclick")
  }
}

#' Plot infection date density for chosen case
#' @param singleRecord Posterior MCMC sample set from TransPhylo
#' @param index Index of case
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param overlayDate Date to overlay on graph
#' @return A vector of posterior infection dates
#' @export
#' @examples
#' plotInfectionDateDensity(record, index, date)
plotInfectionDateDensity <- function(singleRecord, index, burnin=0.5, overlayDate){
  starts <- floor(burnin * length(singleRecord))
  dateList <- lapply(starts:length(singleRecord), function(x) {
    info <- getInfectionDates(singleRecord[[x]]$ctree)[[index]]
    return(info)
  })
  d <- density(as.numeric(dateList))
  plot(d, main=paste0('Infection date density for ', singleRecord[[1]]$ctree$nam[[index]]))
  abline(v=overlayDate, col='red', lwd=2)
}

#' Plot generation time density for chosen case
#' @param singleRecord Posterior MCMC sample set from TransPhylo
#' @param index Index of case
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return A vector of posterior generation times
#' @export
#' @examples
#' plotGenerationTimeDensity(record, index)
plotGenerationTimeDensity <- function(singleRecord, index, burnin=0.5){
  starts <- floor(burnin * length(singleRecord))
  dateList <- lapply(starts:length(singleRecord), function(x) {
    info <- getGenerationTimes(singleRecord[[x]]$ctree)[[index]]
    return(info)
  })
  d <- density(as.numeric(dateList))
  plot(d, main=paste0('Generation time density for ', singleRecord[[1]]$ctree$nam[[index]]))
}

#' Plot time to sampling (from infection) for chosen case
#' @param singleRecord Posterior MCMC sample set from TransPhylo
#' @param index Index of case
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @return A vector of posterior times
#' @export
#' @examples
#' plotSampleTimeDensity(record, index)
plotSampleTimeDensity <- function(singleRecord, index, burnin=0.5){
  starts <- floor(burnin * length(singleRecord))
  dateList <- lapply(starts:length(singleRecord), function(x) {
    info <- getTimesToSampling(singleRecord[[x]]$ctree)[[index]]
    return(info)
  })
  d <- density(as.numeric(dateList))
  plot(d, main=paste0('Sample time density for ', singleRecord[[1]]$ctree$nam[[index]]))
}

#' Plot histogram of generation times for each cluster
#' @param cls_record Posterior MCMC sample set from TransPhylo for all clusters
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param ... parameters to be passed to grid.arrange
#' @return NULL
#' @export
plot_gen_times <- function(cls_record, burnin=0.5, ...){
  gs <- list()
  for(nm in names(cls_record)){
    ends <- length(cls_record[[nm]])
    starts <- floor(burnin * ends)
    gen_times <- map(cls_record[[nm]][starts:ends], ~ getGenerationTimes(.$ctree))
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

#' Plot histogram of times to sampling for each cluster
#' @param cls_record Posterior MCMC sample set from TransPhylo for all clusters
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param ... parameters to be passed to grid.arrange
#' @return NULL
#' @export
plot_times_to_samp <- function(cls_record, burnin=0.5, ...){
  gs <- list()
  for(nm in names(cls_record)){
    ends <- length(cls_record[[nm]])
    starts <- floor(burnin * ends)
    times_to_samp <- map(cls_record[[nm]][starts:ends], ~ getTimesToSampling(.$ctree))
    gs[[nm]] <- ggplot(data.frame(times_to_samp = unlist(times_to_samp)), aes(times_to_samp)) +
      geom_histogram(bins = 20) + ggtitle(nm)
  }
  gridExtra::grid.arrange(grobs = gs, ...)
}

#' Plot histogram of the number of unsampled cases for each cluster
#' @param cls_record Posterior MCMC sample set from TransPhylo for all clusters
#' @param burnin Proportion of the MCMC output to be discarded as burnin
#' @param ... parameters to be passed to grid.arrange
#' @return NULL
#' @export
plot_unsampled_cases <- function(cls_record, burnin=0.5, ...){
  gs <- list()
  for(nm in names(cls_record)){
    ends <- length(cls_record[[nm]])
    starts <- floor(burnin * ends)
    unsampled_cases <- map(cls_record[[nm]][starts:ends], ~ getNumberUnsampled(.$ctree))
    gs[[nm]] <- ggplot(data.frame(unsampled_cases = unlist(unsampled_cases)), aes(unsampled_cases)) +
      geom_histogram(bins = 20) + ggtitle(nm)
  }
  gridExtra::grid.arrange(grobs = gs, ...)
}




