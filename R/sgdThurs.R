#' @export
sgdThurs <- function(data, mu, sigma, rate=0.1, maxiter=1000, tol=1e-9, start, decay=1.1){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0

  #rate is the step size/learning rate

  #helper function
  #return the value of the function to be minimized
  #as well as the gradient w.r.t. whichObs-th observation
  targetThurs <- function(whichObs, scores, data, mu, inv_sigma){
    #let m be the number of varieties,
    #let n be the number of farmers.
    #data is an n*m matrix,
    #data(i, j) represents the rank of variety i by farmer j
    #the entry where varieties are not included is 0

    #(mu, sigma) are parameters of the normal prior

    nFarmers<- nrow(data)
    nVarieties <- ncol(data)
    colnames(data) <- 1:nVarieties #assign labels to varieties
    #initialize
    target_value <- as.numeric(0.5 * (t(scores - mu) %*% inv_sigma %*% (scores - mu)))
    gradient <- 1 / nFarmers * inv_sigma %*% (scores - mu)

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

          #calculate the term involved with cdf of std. normal
          quant <- (scores[win] - scores[lose]) / sqrt(2)
          cdf_term <- pnorm(quant)

          #update the value of the target function
          target_value <- target_value - log(cdf_term)

          grad_change <- (1 / cdf_term) * (1 / sqrt(2 * pi)) *
            exp(-0.5 * quant^2) * (1/ sqrt(2))
          gradient[win] <- gradient[win] - grad_change / nFarmers
          gradient[lose] <- gradient[lose] + grad_change / nFarmers
        }
      }
    }
    return(list(value=target_value, gradient=gradient))
  }#END targetThurs


  #initialize
  colnames(data) <- 1:nVarieties #assign labels to varieties
  start <- scores <- rnorm(nVarieties)
  inv_sigma <- solve(sigma)

  #initialize
  niter <- 0
  nTargetWorse <- 0
  targets <- NULL
  gradients <- NULL
  flag <- TRUE
  #loop until the convergence criteria are met
  while(flag){
    niter <- niter + 1
    res_temp <- targetThurs(i, scores, data, mu, inv_sigma)
    targets <- c(targets, res_temp$value)
    scores <- scores - rate * res_temp$gradient
    gradients <- rbind(gradients, res_temp$gradient)

    #check the convergence criteria: square of the change of target values
    if(niter > 1){
      if((targets[niter] - targets[niter - 1]) > 0){
        nTargetWorse <- nTargetWorse + 1
        rate <- rate / decay
      }
      if((targets[niter] - targets[niter - 1]) ^ 2 < tol | niter > maxiter){
        flag <- FALSE
      }
    }#END niter > 1
  }#END while flag

  return(list(value=targets[niter], niter=niter, scores=scores, startScores=start, targets=targets, gradients=gradients, nTargetWorse=nTargetWorse))
}#END sgdThurs