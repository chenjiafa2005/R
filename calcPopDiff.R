### Function from ploysat package v1.5 ######
calcPopDiff <- function ( freqs, metric, pops = row.names(freqs), 
                          loci = unique(as.matrix(as.data.frame(strsplit(colnames(freqs), split = ".", fixed = TRUE),stringsAsFactors = FALSE))[1,]) ){
  if (!metric %in% c("Fst", "Gst", "Jost's D")) {
    stop("metric must be Fst, Gst, or Jost's D")
  }
  loci <- loci[loci != "Genomes"]
  result <- matrix(0, nrow = length(pops), ncol = length(pops), 
                   dimnames = list(pops, pops))
  if ("Genomes" %in% names(freqs)) {
    genomes <- freqs$Genomes
    names(genomes) <- pops
    GbyL <- FALSE
  } else {
    GbyL <- TRUE
  }
  for (m in 1:length(pops)) {
    for (n in m:length(pops)) {
      hets <- array(0, dim = c(length(loci), 4), dimnames = list(loci, 
                    c("HT", "HS", "HTest", "HSest")))
      if (!GbyL) {
        genomesM <- genomes[pops[m]]
        genomesN <- genomes[pops[n]]
      }
      for (L in loci) {
        if (GbyL) {
          genomesM <- freqs[pops[m], paste(L, "Genomes", sep = ".")]
          genomesN <- freqs[pops[n], paste(L, "Genomes", sep = ".")]
        }
        thesefreqs <- freqs[c(pops[m], pops[n]), grep(paste(L,".", sep = ""), names(freqs), fixed = TRUE)]
        thesefreqs <- thesefreqs[, names(thesefreqs) != paste(L, "Genomes", sep = ".")]
        if (metric == "Fst") {
          avgfreq <- (thesefreqs[1,] * genomesM + thesefreqs[2,] * genomesN)/(genomesM + genomesN)
          hets[L, "HS"] <- ((1 - sum(thesefreqs[1, ]^2)) * genomesM + (1 - sum(thesefreqs[2, ]^2)) * genomesN)/(genomesM + genomesN)
        }
        if (metric %in% c("Jost's D", "Gst")) {
          avgfreq <- (thesefreqs[1, ] + thesefreqs[2,])/2
          hets[L, "HS"] <- (1 - sum(thesefreqs[1, ]^2) + 1 - sum(thesefreqs[2, ]^2))/2
        }
        hets[L, "HT"] <- 1 - sum(avgfreq^2)
        if (metric %in% c("Jost's D", "Gst")) {
          meanGenomes <- 2/(1/genomesM + 1/genomesN)
          hets[L, "HSest"] <- hets[L, "HS"] * meanGenomes/(meanGenomes - 1)
          hets[L, "HTest"] <- hets[L, "HT"] + hets[L, "HSest"]/(2 * meanGenomes)
        }
      }
      if (metric == "Fst") {
        HT <- mean(hets[, "HT"], na.rm=TRUE)  #Jiafa: treat missing value
        HS <- mean(hets[, "HS"], na.rm=TRUE)  #Jiafa: treat missing value
        result[m, n] <- result[n, m] <- (HT - HS)/HT
      }
      if (metric == "Gst") {
        G <- (hets[, "HTest"] - hets[, "HSest"])/hets[,"HTest"]
        result[m, n] <- result[n, m] <- mean(G, na.rm=TRUE)  #Jiafa: treat missing value
      }
      if (metric == "Jost's D") {
        D <- 2 * (hets[, "HTest"] - hets[, "HSest"])/(1 - hets[, "HSest"])
        result[m, n] <- result[n, m] <- mean(D, na.rm=TRUE)  #Jiafa: treat missing value
      }
    }
  }
  return(result)
}
##############################################################################
#
