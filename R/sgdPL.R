#' @export
sgdPL <- function(data, sigma=diag(ncol(data)), rate=0.01, maxiter=1000, tol=1e-8, startVar=1, startScores=scale(rnorm(ncol(data))), decay=1.1){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #for an observation where a variety was not evaluated, its rank is 0

  #(0, sigma) are parameters of the normal prior

  #rate is the step size/learning rate

  #helper function
  #return the value of the function to be minimized
  #the first value of the scores is the variance of the scores
  targetPL <- function(scores, data, inv_sigma){
    sig2 <- scores[1]
    scores <- scores[-1]
    nFarmers <- nrow(data)
    nVarieties <- ncol(data)
    colnames(data) <- 1:nVarieties #assign labels to varieties
    #initialize
    sTSigInvs <- 0.5 * c(t(scores) %*% inv_sigma %*% scores / sig2)
    target_value <- 0.5 * nVarieties * log(sig2) + sTSigInvs #Warning: likelihood sign reversed (trying to minimize...)
    gradient <- inv_sigma %*% scores / sig2 #Warning: gradient sign reversed (will subtract)
    gradSig <- (0.5 * nVarieties - sTSigInvs) / sig2

    #loop over all observations
    for(i in 1:nFarmers){
      #calculate the ranking of the form A>B>C...
      ranks <- data[i, ][data[i, ] != 0]
      ranking <- as.numeric(names(sort(ranks)))

      #the length of i-th observation
      nrank <- length(ranks)

      #loop over all pairwise comparisons
      for(j in 1:(nrank - 1)){
        sum_temp <- 0
        for(k in (j + 1):nrank){
          #ranking[j] wins over ranking[k]
          win <- ranking[j]
          lose <- ranking[k]

          sum_temp <- sum_temp + exp(scores[lose] - scores[win])
          #update the value of the target function
          target_value <- target_value + log(1 + sum_temp)

          after_j <- ranking[(j + 1):nrank] #varieties ranked after variety j
          gradient[after_j] <- gradient[after_j] + (1 / (1 + sum_temp)) * exp(scores[after_j] - scores[win])
          gradient[win] <- gradient[win] - (1/ (1 + sum_temp)) * sum_temp
        }
      }
    }
    gradient <- gradient - mean(gradient) #ensure that the scores always center on 0
    return(list(value=target_value, gradient=c(gradSig, gradient)))
  }#END targetPL

  #initialize
  maxRate <- 2 * rate
  scores <- c(startVar, startScores) #also estimating the variance among scores
  inv_sigma <- solve(sigma)

  #initialize
  niter <- 0
  nTargetWorse <- 0
  targets <- NULL
  parmVals <- scores
  gradients <- NULL
  rates <- NULL
  flag <- TRUE
  #loop until the convergence criteria are met
  while(flag){
    niter <- niter + 1
    res_temp <- targetPL(scores, data, inv_sigma)
    targets <- c(targets, res_temp$value)
    scoreVar <- scores[1]
    scores <- scores - rate * res_temp$gradient
    #prevent dramatic changes in the variance
    varChangeRatio <- scores[1] / scoreVar
    if (varChangeRatio < 0.9) scores[1] <- scoreVar * (1.0 - runif(1, 0, 0.1))
    if (varChangeRatio > 1.1) scores[1] <- scoreVar * (1.0 + runif(1, 0, 0.1))
    parmVals <- rbind(parmVals, scores)
    gradients <- rbind(gradients, res_temp$gradient)
    rates <- c(rates, rate)
    #check the convergence criteria
    if(niter > 1){
      if(targets[niter] > targets[niter - 1]){
        nTargetWorse <- nTargetWorse + 1
        rate <- rate / decay
      } else rate <- min(maxRate, rate * decay) #allow the rate to go up but not too far
      if((targets[niter] - targets[niter - 1]) ^ 2 < tol | niter > maxiter){
        flag <- FALSE
      }
    }#END niter > 1
  }#END while flag
  return(list(value=targets[niter], niter=niter, scoreVar=scores[1], scores=scores[-1], targets=targets, parmVals=parmVals, gradients=gradients, nTargetWorse=nTargetWorse, rates=rates))
}#END sgdPL