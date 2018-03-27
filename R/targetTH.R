#' Likelihood and gradient for the Thurstone model
#' @param scores the scores of the varieties. score[1] is the variance among scores
#' @param data a n * m matrix,
#' where n is the number of observers and m is the number of items to rank;
#' each row vector is a partial ranking (triple comparisons),
#' with i-th element being the rank assigned to item i;
#' if an item is not ranked its partial ranking is 0
#' @param inv_sigma the inverse of the score covariance matrix among varieties
targetTH <- function(scores, data, inv_sigma){
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
    ranks <- data[i, ][data[i, ] != 0]
    ranking <- as.numeric(names(sort(ranks)))

    #the length of i-th observation
    nrank <- length(ranks)

    #loop over all pairwise comparisons
    for(j in 1:(nrank - 1)){
      for(k in (j + 1):nrank){
        #ranking[j] wins over ranking[k]
        win <- ranking[j]
        lose <- ranking[k]

        #calculate the term involved with cdf of std. normal
        quant <- (scores[win] - scores[lose]) / sqrt(2)
        cdf_term <- pnorm(quant)
        #update the value of the target function
        target <- target + log(cdf_term)

        grad_change <- (1 / cdf_term) * (1 / sqrt(2 * pi)) * exp(-0.5 * quant^2) * (1 / sqrt(2))
        gradient[win] <- gradient[win] + grad_change
        gradient[lose] <- gradient[lose] - grad_change
      }
    }
  }
  gradient <- gradient - mean(gradient) #ensure that the scores always center on 0
  return(list(target=target, gradient=c(gradSig, gradient)))
}#END targetTH
