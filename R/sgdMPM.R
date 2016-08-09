#' @export
sgdMPM = function(data, mu, sigma, rate, maxiter = 1000, tol = 1e-9, start, decay = 1.1){

  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0

  #(mu, sigma) is the parameters of the normal prior

  #rate is the starting step size/learning rate
  #start is the starting value of the parameter, ordered as
  #c(score, uncertainty, adherence)


  #helper function

  #transform the ranking data into a count matrix
  rank2count = function(data){

    #inpute data is a vector,
    #where i-th element in the vector is the rank assigned to the i-th item
    #the entry where item not ranked is replaced by 0

    #assume there are no ties

    nitem = length(data)
    names(data) = 1:length(data) #assign labels

    rank = data[data != 0]
    nitem_ranked = length(rank)
    C = matrix(0, nitem, nitem)

    #loop over all pairwise comparisons
    for(i in 1:(nitem_ranked - 1)){

      for(j in (i + 1):nitem_ranked){

        entry = rank[i] - rank[j]
        item_i = as.numeric(names(rank)[i])
        item_j = as.numeric(names(rank)[j])

        if(entry > 0){ #item_j wins over item_i
          C[item_j, item_i] = entry
        } else{  #item_i wins over item_j
          C[item_i, item_j] = -entry
        }

      }
    }

    return(C)

  }








  #helper function

  #return the value of the function to be minimized
  #as well as the gradient w.r.t. index-th observation
  targetMPM = function(index, score, uncertainty, adherence, data, mu, sigma){


    #let m be the number of varieties,
    #let n be the number of farmers.
    #data is an n*m matrix,
    #data(i, j) represents the rank of variety i by farmer j
    #the entry where varieties are not included is 0

    #(mu, sigma) are parameters of the normal prior

    nobs = nrow(data)
    nvar = ncol(data)
    colnames(data) = 1:nvar #assign labels to varieties

    #the first nvar elements are the gradient for score,
    #the next nvar elements are the gradient for uncertainty,
    #the last nobs elements are the gradient for adherence
    gradient = rep(0, (nvar + nvar + nobs))

    #initialize
    inv_sigma = solve(sigma)
    target = as.numeric(0.5 * (t(score - mu) %*% inv_sigma %*% (score - mu)))
    gradient[1:nvar] = 1 / nobs * inv_sigma %*% (score - mu)
    ranks = data[index, ][data[index, ]!= 0] #the index-th observation
    ranking = as.numeric(names(sort(ranks))) #the index-th ranking
    nrank = length(ranking)

    #pre-compute numbers (pairwise) for later computations
    pair_score = outer(score, score, '-')
    pair_uncer = outer(uncertainty, uncertainty, '+')
    frac_term = pair_score / pair_uncer

    ##calculate the target value
    vec1 = rep(0, nobs)
    sum_count = rep(0, nobs)
    for(i in 1:nobs){ #loop over all observations
      count = rank2count(data[i, ]) #count matrix
      sum_count[i] = sum(count)

      temp_res = adherence[i] * frac_term

      #update target value
      target = target - sum(count * temp_res)

      temp_res = exp(temp_res)
      diag(temp_res) = 0

      #store a specific temp_res for later computations
      if(i == index){
        exp_term = temp_res
        mycount = count
      }

      vec1[i] = sum(temp_res)
    }
    #update target value
    target = target + sum(sum_count * log(vec1))


    ##compute gradient
    vec2 = rep(0, nvar)
    vec3 = rep(0, nvar)
    for(i in 1:nvar){
      temp_res = (exp_term[i, ] - exp_term[, i]) * adherence[index]
      vec2[i] = sum(temp_res / pair_uncer[i, ])
      vec3[i] = sum(temp_res * pair_score[i, ] / (pair_uncer[i, ])^2)

      #update gradient w.r.t. score
      gradient[i] = gradient[i] + sum(mycount[, i] * adherence[index] / pair_uncer[i, ])
      gradient[i] = gradient[i] - sum(mycount[i, ] * adherence[index] / pair_uncer[i, ])

      #update gradient w.r.t. uncertainty
      gradient[nvar + i] = gradient[nvar + i] + sum(mycount[, i]
                                                    * adherence[index] * pair_score[, i]
                                                    / (pair_uncer[i, ])^2)
      gradient[nvar + i] = gradient[nvar + i] + sum(mycount[i, ]
                                                    * adherence[index] * pair_score[i, ]
                                                    / (pair_uncer[i, ])^2)

    }
    #update gradient w.r.t. score
    gradient[1:nvar] = gradient[1:nvar] + (1 / vec1[index]) * vec2 * sum_count[index]

    #update gradient w.r.t. uncertainty
    gradient[(nvar + 1):(2 * nvar)] = gradient[(nvar + 1):(2 * nvar)] -
      (1 / vec1[index]) * vec3 * sum_count[index]

    #update gradient w.r.t. adherence
    gradient[2 * nvar + index] = gradient[2 * nvar + index] - sum(mycount * frac_term)
    temp_num = sum(exp_term * frac_term)
    gradient[2 * nvar + index] = gradient[2 * nvar + index] +
      (1 / vec1[index]) * temp_num * sum_count[index]

    return(list(value = target, gradient = gradient))

  }




  #main function starts here
  nobs = nrow(data)
  nvar = ncol(data)
  colnames(data) = 1:nvar #assign labels to varieties
  inv_sigma = solve(sigma)

  #initialize
  niter = 0
  #the first nvar element is the score
  #the next nvar element is the uncertainty
  #the last nobs element is the adherence
  param = start
  target = rep(0, niter)

  flag = TRUE
  #loop until the convergence criteria are met
  while(flag){

    for(i in 1:nobs){

      niter = niter + 1
      score_temp = param[1:nvar]
      uncertainty_temp = param[(nvar + 1):(2*nvar)]
      adherence_temp = param[(2*nvar + 1):(2*nvar + nobs)]
      #evaluate the log-posterior as well as the gradient
      #only used for small dataset (where we want to decide the learning rate)
      #if used for big dataset, where we don't want to
      #evaluate log-posterior everytime, the function should be modified
      res_temp = targetMPM(i, score_temp, uncertainty_temp, adherence_temp, data, mu, sigma)

      #store the value of the target function
      target[niter] = res_temp[[1]]

      #extract the gradient
      gradient = res_temp[[2]]

      #update the parameters
      param = param - rate * gradient

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

  }

  return(list(value = target, niter = niter, score = param[1:nvar],
              uncertainty = param[(nvar + 1):(2 * nvar)],
              adherence = param[(2 * nvar + 1):(2 * nvar + nobs)]))

}
