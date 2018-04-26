#' Gradient descent function generalizable across ranking models
#' @param data a n * m matrix,
#' where n is the number of observers and m is the number of items to rank;
#' each row vector is a partial ranking (triple comparisons),
#' with i-th element being the rank assigned to item i;
#' the entry where that item is not ranked in the partial ranking is replaced by 0
#' @param targetFN calculates the likelihood of the current score and the gradient
#' to move to the next iteration
#' @param sigma the covariance matrix among varieties
#' @param rate the learning rate
#' @param maxiter the maximum number of iterations
#' @param tol the tolerance
#' @param startVar initial guess at variance of scores
#' @param startScores initial guesses at scores
#' @param decay how fast the learning rate decays when likelihood doesn't increase
gradientDescent <- function(data, targetFN=targetPL, sigma=diag(ncol(data)), rate=0.01, maxiter=1000, tol=1e-08, startVar=1, startScores=NULL, decay=0.9){
  on.exit(saveRDS(list(target=targets[nIter], nIter=nIter, scoreVar=scoreVar, scores=scores, startVar=startVar, startScores=startScores, targets=targets, allScores=allScores, allGradients=allGradients, nTargetWorse=nTargetWorse, rates=rates), file=paste("gradDescOut", as.integer(Sys.time()),".rds", sep="")) # For forensics
  # m is the number of varieties,
  # n is the number of farmers.
  # data is an n*m matrix,
  # data(i, j) represents the rank of farmer i by variety j
  # for an observation where a variety was not evaluated, its rank is 0
  # (0, sigma) are parameters of the normal prior
  #rate is the step size = learning rate
  nIter <- 1
  nTargetWorse <- 0
  rates <- NULL
  # Start with reasonable values for the scores
  if (is.null(startScores)){
    startScores <- -scale(colMeans(data, na.rm=T))
  }
  scores <- c(startVar, startScores)
  allScores <- matrix(scores, nrow=1)
  inv_sigma <- solve(sigma)
  # Prepare for the first iteration
  res_temp <- targetFN(scores, data, inv_sigma)
  targets <- res_temp$target
  allGradients <- matrix(res_temp$gradient, nrow=1)
  scores <- scores + rate * res_temp$gradient
  keepGoing <- TRUE
  while(keepGoing){
    res_temp <- targetFN(scores, data, inv_sigma)
    if(res_temp$target < targets[nIter]){ # The likelihood got worse
      # Decrease the rate and revert back to earlier state
      nTargetWorse <- nTargetWorse + 1
      rate <- rate * decay
      scores <- allScores[nrow(allScores),]
    } else{ # The likelihood got better
      # Store the better scores and their gradients
      targets <- c(targets, res_temp$target)
      nIter <- nIter + 1
      allScores <- rbind(allScores, scores)
      allGradients <- rbind(allGradients, res_temp$gradient)
      # Keep going if the likelihood got measurably better
      keepGoing <- (targets[nIter] - targets[nIter - 1]) ^ 2 > tol & nIter < maxiter
    }
    scoreVar <- scores[1]
    scores <- scores + rate * allGradients[nrow(allGradients),]
    # Prevent dramatic changes in the variance. Empirically determined limits
    varChangeRatio <- scores[1] / scoreVar
    if (varChangeRatio < 0.6) scores[1] <- scoreVar * 0.6
    if (varChangeRatio > 1.4) scores[1] <- scoreVar * 1.4
    if (scores[1] < 0.01) scores[1] <- 0.01
    if (scores[1] > 10) scores[1] <- 10
    rates <- c(rates, rate)
  }#END while keepGoing
  # NOTE all the ranking methods got ranking flipped, so flip it back
  scoreVar <- scores[1]
  scores <- -scores[-1]
  # on.exit()
  return(list(target=targets[nIter], nIter=nIter, scoreVar=scoreVar, scores=scores, startVar=startVar, startScores=startScores, targets=targets, allScores=allScores, allGradients=allGradients, nTargetWorse=nTargetWorse, rates=rates))
}#END gradientDescent
