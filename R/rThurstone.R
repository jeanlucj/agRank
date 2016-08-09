#' @export
rThurstone = function(S, Svar){
  #S: Score vector
  #Svar: variance vector
  nitem = length(S)
  utility = rnorm(nitem, S, sqrt(Svar))
  names(utility) = 1:nitem #assign labels to varieties
  #ranking is the items listed in the order A succeeds B succeeds C, etc.
  ranking = as.numeric(names(sort(utility, decreasing = TRUE)))

  #calculate ranks of each item
  ranks = match(1:nitem, ranking)

  return(list(ranks = ranks, ranking = ranking))
}
