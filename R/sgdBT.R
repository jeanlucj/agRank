#' @export
sgdBT = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start, decay = 1.1){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0

  #rate is the step size/learning rate

  #helper function
  #return the value of the function to be minimized
  #as well as the gradient w.r.t. index-th observation
  targetBT = function(index, score, adherence, data, mu, sigma){

    #let m be the number of varieties,
    #let n be the number of farmers.
    #data is an n*m matrix,
    #data(i, j) represents the rank of variety i by farmer j
    #the entry where varieties are not included is 0

    #(mu, sigma) are parameters of the normal prior

    nobs = nrow(data)
    nvar = ncol(data)
    colnames(data) = 1:nvar #assign labels to varieties

    #the first nvar element is the gradient for score
    #the last nobs element is the gradient for adherence
    gradient = rep(0, nvar + nobs)

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

          exp_term = exp(-adherence[i] * (score[win] - score[lose]))

          #update the value of the target function
          target_value = target_value + log(1 + exp_term)

          if(index == i){
            #update the gradient w.r.t. score
            gradient[win] = gradient[win] - adherence[i] * exp_term / (1 + exp_term)
            gradient[lose] = gradient[lose] + adherence[i] * exp_term / (1 + exp_term)

            #update the gradient w.r.t. adherence
            gradient[nvar + i] = gradient[nvar + i] +
              (-score[win] + score[lose]) * exp_term / (1 + exp_term)

          }

        }
      }


    }

    return(list(value = target_value, gradient = gradient))
  }






  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  inv_sigma = solve(sigma)

  #initialize
  #the first nvar element is the score
  #the last nobs element is the adherence
  param = start
 
ADAGRAD = function(data, maxiter=10000){
  #convert data.frame dataSet in matrix
  data = matrix(unlist(data), ncol=ncol(data), byrow=FALSE)
  #parse data into input and output
  inputData = data[,1:ncol(data)-1]
  outputData = data[,ncol(data)]
  
  #initiaize matrices
  temporaryparam = matrix(0,ncol=length(param), nrow=1)
  updateRule = matrix(0, ncol=length(param), nrow=1)
  gradientList = matrix(NA,nrow=1, ncol=0)
  stochasticList = sample(1:nrow(data), maxiter, replace=TRUE)
  target=c(0,maxiter)
  #loop the gradient descent
  for(niter in 1:maxiter){
    error = (inputData[stochasticList[niter],] %*% (param)) - outputData[stochasticList[niter]]
    for(column in 1:length(param)){
      #calculate gradient
      res_temp = targetBT(column, score_temp, adherence_temp, data, mu, sigma)
      target[niter] = res_temp[[1]]
      gradient = as.numeric(error) * t(as.matrix(res_temp[[2]]))
      #adagrad update rule calculation
      gradientList = cbind(gradientList, gradient)
      gradientSum = as.numeric(sqrt(gradientList %*% t(gradientList)))
      updateRule = as.numeric(rate / gradientSum) * gradient
      temporaryparam = param - rate*updateRule
       #check the convergence criteria: square of the change of target values
       if(niter > 1){
        if((target[niter] - target[niter - 1]) ^ 2 < tol | niter > maxiter){
          flag = FALSE
          break
        }

        #update learning rate if the target value don't decrease
        if((target[niter - 1] - target[niter]) / target[niter - 1] < 0){
          rate = rate / decay
        }
      }
    }
    #update all theta in the current iteration
    param = temporaryparam
  }
  
  return(list(value = target, niter = niter, score = param[1:nvar],
              adherence = param[(nvar + 1):(nvar + nobs)]))
  }
  return(list(value = target, niter = niter, score = param[1:nvar],
              adherence = param[(nvar + 1):(nvar + nobs)]))
}
