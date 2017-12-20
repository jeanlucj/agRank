#' @export
sgdBT <- function(data, sigma=diag(ncol(data)), rate=0.001, maxiter=1000, tol=1e-9, start=rnorm(ncol(data)), decay=1.1){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0

  #(0, sigma) are parameters of the normal prior

  #rate is the step size/learning rate

  #helper function
  #return the value of the function to be minimized
  #CH the first value of the scores is the variance
  targetBT <- function(scores, data, inv_sigma){
    sig2 <- scores[1] #CH
    scores <- scores[-1] #CH
    nFarmers <- nrow(data)
    nVarieties <- ncol(data)
    colnames(data) <- 1:nVarieties #assign labels to varieties
    #initialize
    target_value <- as.numeric(0.5 * (t(scores) %*% inv_sigma %*% scores)) / sig2 #Warning: likelihood sign reversed (trying to minimize...)
    gradient <- inv_sigma %*% scores #Warning: gradient sign reversed (will subtract)
    gradSig <- (nVarieties - target_value) / sig2 #CH boy: I hope this is right...

    #loop over all observations
    for(i in 1:nFarmers){
      #calculate the ranking of the form A>B>C...
      ranks <- data[i, ][data[i, ] !=0]
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
          target_value <- target_value + log(1 + exp_term)

          gradient[win] <- gradient[win] - exp_term / (1 + exp_term)
          gradient[lose] <- gradient[lose] + exp_term / (1 + exp_term)

        }
      }
    }
    gradient <- gradient - mean(gradient) #CH ensure that the scores always center on 0
    return(list(value=target_value, gradient=c(gradSig, gradient)))
  }#END targetBT

  #initialize
  nFarmers <- nrow(data)
  nVarieties <- ncol(data)
  colnames(data) <- 1:nVarieties #assign labels to varieties
  scores <- scale(rnorm(nVarieties)) #CH also estimating the variance among scores
  start <- scores <- c(1, scores) #CH also estimating the variance among scores
  inv_sigma <- solve(sigma)

  #initialize
  niter <- 0
  nTargetWorse <- 0
  targets <- NULL
  gradients <- NULL
  rates <- NULL
  flag <- TRUE
  #loop until the convergence criteria are met
  while(flag){
    niter <- niter + 1
    res_temp <- targetBT(scores, data, inv_sigma)
    targets <- c(targets, res_temp$value)
    scores <- scores - rate * res_temp$gradient
    gradients <- rbind(gradients, res_temp$gradient)
    rates <- c(rates, rate)
    #check the convergence criteria: square of the change of target values
    if(niter > 1){
      if(targets[niter] > targets[niter - 1]){
        nTargetWorse <- nTargetWorse + 1
        rate <- rate / decay
      } else rate <- min(1, rate * (1 + (decay - 1) / 2)) #CH also allow the rate to go up slowly
      if((targets[niter] - targets[niter - 1]) ^ 2 < tol | niter > maxiter){
        flag <- FALSE
      }
    }#END niter > 1
  }#END while flag
saveRDS(list(value=targets[niter], niter=niter, scoreVar=scores[1], scores=scores[-1], startScores=start, targets=targets, gradients=gradients, nTargetWorse=nTargetWorse, rates=rates), file="sgdBTout.RDS") #CH for forensics
  return(list(value=targets[niter], niter=niter, scoreVar=scores[1], scores=scores[-1], startScores=start, targets=targets, gradients=gradients, nTargetWorse=nTargetWorse))
}#END sgdBT
