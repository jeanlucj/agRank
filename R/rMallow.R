library(stats) #sampler from multinomial
#' @import stats
#' @export
rMallow = function(sigma, phi){
  #sigma is the reference ranking in the form (sig1, sig2,..., sigm )
  #phi is the dispersion parameter
  if(phi <= 0 | phi > 1){
    warning('phi should be in (0, 1]')
  }
  nitem = length(sigma)
  #initialize
  ranking = c(sigma[1])

  #repeated insersion
  for(i in 2:nitem){
    j = 1:i #the place to insert
    if(phi == 1){
      pvec = rep(1/i, i) #insersion probabilities
    } else{
      pvec = phi^(i-j) * (1-phi) / (1-phi^i) #insersion probabilities
    }
    ind = which(rmultinom(1, 1, pvec) == 1) - 1 #insert after which index
    ranking = append(ranking, sigma[i], ind)
  }

  #calculate ranks of each item
  ranks = match(1:nitem, ranking)

  return(list(ranks = ranks, ranking = ranking))

}

