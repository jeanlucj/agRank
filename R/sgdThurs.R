#' @export
sgdThurs = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start, decay){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0

  #rate is the step size/learning rate

  #helper function
  #return the value of the function to be minimized
  #as well as the gradient w.r.t. index-th observation
  targetThurs = function(index, score, adherence, data, mu, sigma){
    #let m be the number of varieties,
    #let n be the number of farmers.
    #data is an n*m matrix,
    #data(i, j) represents the rank of variety i by farmer j
    #the entry where varieties are not included is 0

    #(mu, sigma) are parameters of the normal prior

    nobs = nrow(data)
    nvar = ncol(data)
    

    #the first nvar element is the gradient for score
    #the last nobs element is the gradient for adherence
    gradient = rep(0, nvar)
    #initialize
    inv_sigma = solve(sigma)
    target_value = as.numeric(0.5 * (t(score - mu) %*% inv_sigma %*% (score - mu)))
    gradient[1:nvar] = 1 / nobs * inv_sigma %*% (score - mu)

    #loop over all observations
    for(i in 1:nobs){

      #calculate the ranking of the form A>B>C...
      ranks = data[i, ][data[i, ] != 0]
      ranking = as.numeric(names(sort(ranks)))

      #the length of i-th observation
      nrank = length(ranks)

      #loop over all pairwise comparisons
      for(j in 1:(nrank - 1)){

        for(k in (j + 1):nrank){

          #ranking[j] wins over ranking[k]
          win = ranking[j]
          lose = ranking[k]

          #calculate the term involved with cdf of std. normal
          quant = sqrt(adherence[i]) * (score[win] - score[lose]) / sqrt(2)
          cdf_term = pnorm(quant)

          #update the value of the target function
          target_value = target_value - log(cdf_term)

          if(i == index){
            #update the gradient w.r.t. score
            grad_change = (1 / cdf_term) * (1 / sqrt(2 * pi)) *
              exp(-0.5 * quant^2) * (sqrt(adherence[i]) / sqrt(2))
            gradient[win] = gradient[win] - grad_change
            gradient[lose] = gradient[lose] + grad_change

            #update the gradient w.r.t. adherence
            gradient[nvar + i] = gradient[nvar + i] -
              (1 / cdf_term) * (1 / sqrt(2 * pi)) * exp(-0.5 * quant^2) *
              (score[win] - score[lose]) / (sqrt(2) * 2 * sqrt(adherence[i]))

          }


        }
      }
}
    }

    return(gradient)

  }






  param = start
 
ADAGRAD = function(data, alpha=0.1, maxiter=10){
  #convert data.frame dataSet in matrix
  dataTrain = matrix(unlist(data), ncol=ncol(data), byrow=FALSE)
  #initialize theta
  param= rep(1,ncol(data))
  #bind 1 column to dataTrain
  dataTrain <- cbind(1, dataTrain)
  inputData <- dataTrain[,1:ncol(dataTrain)-1]
  outputData <- dataTrain[,ncol(dataTrain)]
  #temporary variables
  temporaryparam <- matrix(0,ncol=length(param), nrow=1)
  updateRule <- matrix(0, ncol=length(param), nrow=1)
  gradientList <- matrix(NA,nrow=1, ncol=0)
  stochasticList <- sample(1:nrow(dataTrain), 10, replace=TRUE)

  #loop the gradient descent
  for(niter in 1:maxiter){
    for(column in 1:length(param)){
      #calculate gradient
      gradient <- targetThurs(column, score, adherence, data, mu, sigma)[column]
      #adagrad update rule calculation
      gradientList <- cbind(gradientList, gradient)
      gradientSum <- as.numeric(sqrt(gradientList %*% t(gradientList)))
      updateRule[1,column] <- (alpha / gradientSum) * gradient
      temporaryparam[1,column] = param[column] - updateRule[1,column]
    }
    #update all theta in the current iteration
    param <- temporaryparam
  }
  
  return(param)
}
