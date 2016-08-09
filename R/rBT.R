#' @export
rBT = function(S){
  #S: score vector
  nitem = length(S)
  #generate count matrix C
  C = matrix(0, nitem, nitem)

  for(i in 1:(nitem - 1)){
    for(j in (i + 1):nitem){
      #generate indicator of whether v_i defeats v_j
      count = rbinom(1, 1, 1 / (1 + exp(-S[i] + S[j])))
      if(count == 1){
        C[i, j] = 1
      } else{
        C[j, i] = 1
      }
    }
  }

  #number of wins
  wins = apply(C, 1, sum)

  #assign labels to farmers
  names(wins) = 1:nitem

  ranking = as.numeric(names(sort(wins, decreasing = TRUE)))

  #calculate ranks of each item
  ranks = match(1:nitem, ranking)

  return(list(ranks = ranks, ranking = ranking))
}

