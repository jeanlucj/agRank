#' Rank aggregation using linear models
#'
#' Aggregate the ranking from triple comparisons into a consensus ranking using linear models.
#'
#' @param data a n * m matrix,
#' where n is the number of observers and m is the number of items each observer ranked;
#' each row vector is a partial or complete ranking,
#' with i-th element being the rank assigned to item i;
#' the entry where that item is not ranked in the partial ranking is replaced by 0
#' @param K the additive relationship matrix; if provided, it will be used to specify the covariance structure
#' of the linear model
#'
#' @return Return a list with two components:
#'     \item{ranks}{a vector where the i-th element is the rank assigned to the i-th item.}
#'     \item{scores}{a vector of the estimated scores of the varieties.}
#'     \item{scoreVar}{a scalar of the estimated score variance. This variance is method specific.}
#'
#' @export

###this function can only be used to analyze tricot comparisons
#require: each variety at least be compared once, i.e. nobs * 3 >= nvar
rankLM = function(data, K=NA){
  library(lme4)
  library(EMMREML)
  #let m be the number of varieties,
  #let n be the number of farmers.
  #data is an n*m matrix,
  #data(i, j) represents the rank of variety i by farmer j
  #the entry where varieties are not included is 0

  #if K is provided, then the analysis is done assuming it is an additive relationship matrix
  #otherwise, variety scores are assumed to be iid

  data_linear <- which(data > 0, arr.ind = T)
  data_linear <- cbind(data_linear, data[data_linear])
  data_linear <- as.data.frame(data_linear)
  colnames(data_linear) <- c('farmer', 'variety', 'rank')
  data_linear$variety <- factor(data_linear$variety)
  if(!is.matrix(K)){
    fit <- lmer(rank ~ 1 + (1 | variety), data = data_linear)
    scores <- coef(fit)$variety #larger value means less competitive
    scoreVar <- summary(fit)$varcor[1]
  } else{
    y <- data_linear$rank #the response
    X <- as.matrix(rep(1, length(y))) #the intercept
    #Z: design matrix for random effects
    Z <- model.matrix(as.formula("rank ~ -1 + variety"), data=data_linear)
    fit_m <- emmreml(y, X, Z, K)
    scores <- fit_m$uhat
    scoreVar <- fit_m$Vu
  }
  ranks <- rank(scores)
  return(list(ranks = ranks, scores=scores, scoreVar=scoreVar))
}
