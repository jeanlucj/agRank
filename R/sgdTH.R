#' Gradient descent using the Thurstone model
#' @param data a n * m matrix,
#' where n is the number of observers and m is the number of items to rank;
#' each row vector is a partial ranking (triple comparisons),
#' with i-th element being the rank assigned to item i;
#' the entry where that item is not ranked in the partial ranking is replaced by 0
#' @param sigma the covariance matrix among varieties
#' @param rate the learning rate
#' @param maxiter the maximum number of iterations
#' @param tol the tolerance
#' @param startVar initial guess at variance of scores
#' @param startScores initial guesses at scores
#' @param decay how fast the learning rate decays when log post doesn't
#' @export
sgdTH <- function(data, sigma=diag(ncol(data)), rate=0.01, maxiter=1000, tol=1e-8, startVar=1, startScores=scale(rnorm(ncol(data))), decay=0.2){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of farmer i by variety j
  #for an observation where a variety was not evaluated, its rank is 0

  #(0, sigma) are parameters of the normal prior

  #rate is the step size/learning rate

  #helper function
  #return the value of the function to be minimized
  #the first value of the scores is the variance of the scores
  targetTH <- function(scores, data, inv_sigma){
    sig2 <- scores[1]
    scores <- scores[-1]
    nFarmers <- nrow(data)
    nVarieties <- ncol(data)
    colnames(data) <- 1:nVarieties #assign labels to varieties
    #initialize
    sTSigInvs <- -0.5 * c(t(scores) %*% inv_sigma %*% scores) / sig2
    target_value <- sTSigInvs - 0.5 * nVarieties * log(sig2)
    gradient <- -inv_sigma %*% scores / sig2
    gradSig <- (-0.5 * nVarieties - sTSigInvs) / sig2

    #loop over all observations
    for(i in 1:nFarmers){
      #calculate the ranking of the form A>B>C...
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
          target_value <- target_value + log(cdf_term)

          grad_change <- (1 / cdf_term) * (1 / sqrt(2 * pi)) * exp(-0.5 * quant^2) * (1 / sqrt(2))
          gradient[win] <- gradient[win] + grad_change
          gradient[lose] <- gradient[lose] - grad_change
        }
      }
    }
    gradient <- gradient - mean(gradient) #ensure that the scores always center on 0
    return(list(value=target_value, gradient=c(gradSig, gradient)))
  }#END targetTH

  #initialize
  initialRate <- rate
  #also estimating the variance among scores
  scores <- c(startVar, startScores)
  inv_sigma <- solve(sigma)

  #initialize
  nIter <- 0
  nTargetWorse <- 0
  targets <- NULL
  parmVals <- scores
  gradients <- NULL
  rates <- NULL
  keepGoing <- TRUE
  #loop until the convergence criteria are met
  while(keepGoing){
    nIter <- nIter + 1
    res_temp <- targetTH(scores, data, inv_sigma)
    targets <- c(targets, res_temp$value)
    if(nIter > 2){
      if(targets[nIter] < targets[nIter - 1]){
        nTargetWorse <- nTargetWorse + 1
        rate <- rate * decay
      } else{
        scoreVar <- scores[1]
        scores <- scores + rate * res_temp$gradient
        #prevent dramatic changes in the variance
        varChangeRatio <- scores[1] / scoreVar
        if (varChangeRatio < 0.6) scores[1] <- scoreVar * 0.6
        if (varChangeRatio > 1.4) scores[1] <- scoreVar * 1.4
        if (scores[1] < 0.01) scores[1] <- 0.01
        if (scores[1] > 10) scores[1] <- 10
      }
      #check the convergence criteria
      keepGoing <- (targets[nIter] - targets[nIter - 2]) ^ 2 > tol & nIter < maxiter
    } else{
      scores <- scores + rate * res_temp$gradient
    }#END if nIter > 2
    parmVals <- rbind(parmVals, scores)
    gradients <- rbind(gradients, res_temp$gradient)
    rates <- c(rates, rate)
  }#END while keepGoing
  saveRDS(list(value=targets[nIter], nIter=nIter, scoreVar=scores[1], scores=scores[-1], startVar=startVar, startScores=startScores, targets=targets, parmVals=parmVals, gradients=gradients, nTargetWorse=nTargetWorse, rates=rates), file="sgdTHout.RDS") # for forensics
  return(list(value=targets[nIter], nIter=nIter, scoreVar=scores[1], scores=scores[-1], startVar=startVar, startScores=startScores, targets=targets, parmVals=parmVals, gradients=gradients, nTargetWorse=nTargetWorse, rates=rates))
}#END sgdTH
