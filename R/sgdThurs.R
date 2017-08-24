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
 #assign labels to varieties

    #the first nvar element is the gradient for score
    #the last nobs element is the gradient for adherence
    gradient = rep(0, 6)

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

          quant =(score[win] - score[lose]) / sqrt(2)
cdf_term = pnorm(quant)

#update the value of the target function
target_value = target_value - log(cdf_term)

if(i == index){
  #update the gradient w.r.t. score
  grad_change = (1 / cdf_term) * (1 / sqrt(2 * pi)) *
    exp(-0.5 * quant^2)  / sqrt(2)
  gradient[win] = gradient[win] - grad_change
  gradient[lose] = gradient[lose] + grad_change
  
  #update the gradient w.r.t. adherence
  gradient[nvar + i] = gradient[nvar + i] -
    (1 / cdf_term) * (1 / sqrt(2 * pi)) * exp(-0.5 * quant^2) *
    (score[win] - score[lose]) / (sqrt(2) * 2)
  
}


}
}

}
}
return(list(value = target_value, gradient = gradient))

}






  #initialize
  #the first nvar element is the score
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  inv_sigma = solve(sigma)

  #initialize
  
  #the first nvar element is the score
  #the last nobs element is the adherence
  param = start
  target = rep(0, maxiter)
  
  #convert data.frame dataSet in matrix
  data <- matrix(unlist(data), ncol=ncol(data), byrow=FALSE)
  
  #temporary variables
 
  score=ADAGRAD(data)


