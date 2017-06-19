#' @export
sgdBT = function(data, mu, sigma, rate, maxiter = 100, tol = 1e-9, start, decay = 1.1){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0

  #rate is the step size/learning rate

  #helper function
  #return the value of the function to be minimized
  #as well as the gradient w.r.t. index-th observation
  targetBT = function(index, score,  data, mu, sigma){

    #let m be the number of varieties,
    #let n be the number of farmers.
    #data is an n*m matrix,
    #data(i, j) represents the rank of variety i by farmer j
    #the entry where varieties are not included is 0

    #(mu, sigma) are parameters of the normal prior

    nobs = nrow(data)
    nvar = ncol(data)
     #assign labels to varieties

    #the first nvar element is the gradient for score
    #the last nobs element is the gradient for adherence
    gradient = rep(0, nvar)

    #initialize
    inv_sigma = solve(sigma)
    target_value = as.numeric( 0.5*(t(score - mu) %*% inv_sigma %*% (score - mu)))
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

          exp_term = exp (score[win] - score[lose])

          #update the value of the target function
          target_value = target_value + log(1 + exp_term)

          if(index == i){
            #update the gradient w.r.t. score
            gradient[win] = gradient[win] -  exp_term / (1 + exp_term)
            gradient[lose] = gradient[lose] +  exp_term / (1 + exp_term)

            #update the gradient w.r.t. adherence
            gradient[nvar + i] = gradient[nvar + i] +
              (-score[win] + score[lose]) * exp_term / (1 + exp_term)

          }

        }
      }


    }

    return(list(value = target_value, gradient = gradient))
  }







  #initialize
  #the first nvar element is the score
  param = start
  target = rep(1,maxiter)

  flag = TRUE
  #loop until the convergence criteria are met
  while(flag){
    #bind 1 column to dataTrain
  dataTrain <- cbind(1, data)
  #parse dataTrain into input and output
  inputData <- dataTrain[,1:ncol(dataTrain)-1]
  outputData <- dataTrain[,ncol(dataTrain)]
  #temporary variables
  temporaryparam <- matrix(ncol=length(param), nrow=1)
  updateRule <- matrix(0, ncol=length(param), nrow=1)
  gradientList <-  c(NA)
  #constant variables
  rowLength <- nrow(dataTrain)

  stochasticList <- sample(1:nrow(dataTrain), maxiter, replace=TRUE)
  nobs = nrow(dataTrain)
  nvar = ncol(dataTrain)
  
   for(iteration in 1:maxiter){
    
    for(i in 1:nvar){
      #calculate gradient
      score_temp = param[1:nvar]
      res_temp = targetBT(i, score_temp,  data, mu, sigma)
      gradient <-  res_temp[[2]]
      #adagrad update rule calculation
      gradientList <- append(gradientList, gradient)
      gradientSum <- sqrt(gradientList %*% t(gradientList))
      updateRule[1,i] <- (0.1 / gradientSum) * gradient
      temporaryparam[1,i] = param[1,i] - updateRule[1,i]
    }
    #update all theta in the current iteration
    param <- temporaryparam
  }
      #check the convergence criteria: square of the change of target values
      if(niter > 1){
        if((target[iteration] - target[iteration - 1]) ^ 2 < tol | iteration > maxiter){
          flag = FALSE
          break
        }

        #update learning rate if the target value don't decrease
        if((target[iteration - 1] - target[iteration]) / target[iteration - 1] < 0){
          rate = rate / decay
        }
      }


    }
    return(list(value = target, niter = iteration, score = param[1:nvar]))
  }

 
