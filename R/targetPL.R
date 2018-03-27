#' Likelihood and gradient for the Plackett-Luce model
#' @param scores the scores of the varieties. score[1] is the variance among scores
#' @param data a n * m matrix,
#' where n is the number of observers and m is the number of items to rank;
#' each row vector is a partial ranking (triple comparisons),
#' with i-th element being the rank assigned to item i;
#' if an item is not ranked its partial ranking is NA
#' @param inv_sigma the inverse of the score covariance matrix among varieties
targetPL <- function(scores, data, inv_sigma){
  sig2 <- scores[1]
  scores <- scores[-1]
  nFarmers <- nrow(data)
  nVarieties <- ncol(data)
  colnames(data) <- 1:nVarieties #assign labels to varieties
  #initialize
  sTSigInvs <- -0.5 * c(t(scores) %*% inv_sigma %*% scores) / sig2
  target <- sTSigInvs - 0.5 * nVarieties * log(sig2)
  gradient <- -inv_sigma %*% scores / sig2
  gradSig <- (-0.5 * nVarieties - sTSigInvs) / sig2

  #loop over all observations
  for(i in 1:nFarmers){
    ranks <- data[i, !is.na(data[i,])]
    ranking <- as.numeric(names(sort(ranks)))

    #the length of i-th observation
    nrank <- length(ranks)

    #loop over all winners
    for(j in 1:(nrank - 1)){
      win <- ranking[j]
      lose <- ranking[(j + 1):nrank] #varieties ranked after variety j
      allExp <- exp(scores[j:nrank])
      denom <- sum(allExp)
      gradient[lose] <- gradient[lose] - allExp[-1] / denom
      gradient[win] <- gradient[win] + 1 - allExp[1] / denom
      #update the value of the target function
      target <- target + log(allExp[1] / denom)
    }
  }
  gradient <- gradient - mean(gradient) #ensure that the scores always center on 0
  return(list(target=target, gradient=c(gradSig, gradient)))
}#END targetPL
