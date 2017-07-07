#' @export
nobs = nrow(partialRanks)
nvar = ncol(partialRanks)
sgdThurs = function(data, mu, sigma, rate=0.1, maxiter = 1000, tol = 1e-9, start, decay = 1.1){
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0
  
  #rate is the step size/learning rate
  
  #helper function
  #return the value of the function to be minimized
  #as well as the gradient w.r.t. index-th observation
  targetThurs = function(index, score, data, mu, sigma){
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
          quant = (score[win] - score[lose]) / sqrt(2)
          cdf_term = pnorm(quant)
          
          #update the value of the target function
          target_value = target_value - log(cdf_term)
          
          if(i == index){
            #update the gradient w.r.t. score
          grad_change = (1 / cdf_term) * (1 / sqrt(2 * pi)) *
          exp(-0.5 * quant^2) * (1/ sqrt(2))
          gradient[win] = gradient[win] - grad_change
          gradient[lose] = gradient[lose] + grad_change
            
            #update the gradient w.r.t. adherence
            
          }
          
          
        }
      }
      
    }
    
    return(list(value = target_value, gradient = gradient))
    
  }
  
  

  
  #initialize
  #the first nvar element is the score
  #the last nobs element is the adherence
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nVarieties #assign labels to varieties
  inv_sigma = solve(sigma)
  start = rnorm(ncol(data),10,1)
  #initialize
  niter = 0
  #the first nvar element is the score
  #the last nobs element is the adherence
  data = matrix(unlist(data), ncol=ncol(data), byrow=FALSE)
  data <- cbind(1, data)
  #parse data into input and output
  inputData = data[,1:ncol(data)-1]
  outputData = data[,ncol(data)]
  param = start
  target = rep(0, niter)
  gradient= list()
 
  temporaryparam = matrix(0,ncol=length(param), nrow=1)
  updateRule = matrix(0, ncol=length(param), nrow=1)
  gradientList = matrix(NA,nrow=1, ncol=0)
  stochasticList = sample(1:nrow(data), 1000, replace=TRUE)
  flag = TRUE
  #loop until the convergence criteria are met
  while(flag){
    
    for(i in 1:nobs){
      niter = niter + 1
      #error = (param%*%inputData[stochasticList[niter],]) - outputData[stochasticList[niter]]
      score_temp = param[1:nvar]
      
      #evaluate the log-posterior as well as the gradient
      #only used for small dataset (where we want to decide the learning rate)
      #if used for big dataset, where we don't want to
      #evaluate log-posterior everytime, the function should be modified
      res_temp = targetThurs(i,score_temp, partialRank, mu, sigma)
      #store the value of the target function
      target[niter] = res_temp$value
      
      #extract the gradient
      gradient[[niter]] = res_temp[[2]]
      #update the parameters
      gradientList = append(gradientList,gradient[[niter]])
      gradientSum = as.numeric(sqrt(t(gradientList) %*% (gradientList)))
      updateRule = as.numeric(rate/ gradientSum) * gradient[[niter]]
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
      param=temporaryparam
      
    }
    
  }
  
  return(list(value = target, niter = niter, score = param[1:nvar],gradient))
  
}
 

