#' Gradient descent using the Bradley-Terry model
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
sgdBT <- function(data, sigma=diag(ncol(data)), rate=0.01, maxiter=1000, tol=1e-8, startVar=1, startScores=scale(rnorm(ncol(data))), decay=0.2){
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
  targetBT <- function(scores, data, inv_sigma){
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

          #calculate the term involved with exp in the denominator
          exp_term <- exp(scores[lose] - scores[win])
          #update the value of the target function
          target_value <- target_value - log(1 + exp_term)

          gradient[win] <- gradient[win] + exp_term / (1 + exp_term)
          gradient[lose] <- gradient[lose] - exp_term / (1 + exp_term)
        }
      }
    }
    gradient <- gradient - mean(gradient) #ensure that the scores always center on 0
    return(list(value=target_value, gradient=c(gradSig, gradient)))
  }#END targetBT

  #initialize
  initialRate <- rate
  #also estimating the variance among scores
  scores <- c(startVar, startScores)
  inv_sigma <- solve(sigma)

  #initialize
  nIter <- 0
  nTargetWorse <- 0
  targets <- NULL
  allScores <- scores
  allGradients <- NULL
  rates <- NULL
  keepGoing <- TRUE
  #loop until the convergence criteria are met
  while(keepGoing){
    res_temp <- targetBT(scores, data, inv_sigma)
    if (nIter < 2){ # Assume that the likelihood improves at least in the first two iterations
      nIter <- nIter + 1
      targets <- c(targets, res_temp$value)
      scores <- scores + rate * res_temp$gradient
    } else{ # The likelihood got worse, so decrease the rate and revert back to earlier state
      if(targets[nIter] < targets[nIter - 1]){
        nTargetWorse <- nTargetWorse + 1
        initialRate <- initialRate * decay
        rate <- initialRate
        scores <- allScores[nrow(allScores),]
      } else{ # The likelihood got better, so store the better scores and their gradients
        nIter <- nIter + 1
        targets <- c(targets, res_temp$value)
        allScores <- rbind(allScores, scores)
        allGradients <- rbind(allGradients, res_temp$gradient)
        # If the likelihood got worse, always keep going
        # If the likelihood got better, keep going if it got measurably better
        keepGoing <- (targets[nIter] - targets[nIter - 1]) ^ 2 > tol & nIter < maxiter
      }
      scoreVar <- scores[1]
      scores <- scores + rate * allGradients[nrow(allGradients),]
      #prevent dramatic changes in the variance
      varChangeRatio <- scores[1] / scoreVar
      if (varChangeRatio < 0.6) scores[1] <- scoreVar * 0.6
      if (varChangeRatio > 1.4) scores[1] <- scoreVar * 1.4
      if (scores[1] < 0.01) scores[1] <- 0.01
      if (scores[1] > 10) scores[1] <- 10
    }
    rates <- c(rates, rate)
  }#END while keepGoing
  saveRDS(list(value=targets[nIter], nIter=nIter, scoreVar=scores[1], scores=scores[-1], startVar=startVar, startScores=startScores, targets=targets, allScores=allScores, allGradients=allGradients, nTargetWorse=nTargetWorse, rates=rates), file="sgdBTout.RDS") # for forensics
  return(list(value=targets[nIter], nIter=nIter, scoreVar=scores[1], scores=scores[-1], startVar=startVar, startScores=startScores, targets=targets, allScores=allScores, allGradients=allGradients, nTargetWorse=nTargetWorse, rates=rates))
}#END sgdBT
