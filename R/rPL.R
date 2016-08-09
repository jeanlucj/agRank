library(FAdist) #for sampling from Gumbel distribution
#this sampler is based on the equivalence
#of Thurstonian model and Luce's Choice
#when utilities are drawn independently
#' @import FAdist
#' @export
rPL = function(S, beta = 2){
  #browser()
  #S: score vector
  #beta: arbitrary positive scale parameter
  nitem = length(S)

  #compute the location parameter of Gumbel distribution
  mu = beta * S

  #sample from Gumbell distribution
  utility = rgumbel(nitem, scale = beta, location = mu)
  names(utility) = 1:nitem #assign labels to items
  #get ranking
  ranking = as.numeric(names(sort(utility, decreasing = TRUE)))

  #calculate ranks of each item
  ranks = match(1:nitem, ranking)

  return(list(ranks = ranks, ranking = ranking))
}
