#' Simulation results when 10\% outliers are included
#'
#'We do a simulation to see the performance of models when outliers are included.
#'We randomly choose 10\% farmers and the ranking given by those farmers are completely random permutations.
#'We include LM, TH, BT and PL with genotypic information.
#'We run the procedure 100 times with h^2 = 0.2, 0.5, 0.8 with nvar = 240 and nobs = 180, 240, 300, 360.
#'We use the alpha design to allocate varieties to farmers.
#'
#' @docType data
#'
#' @usage data(resOutlier)
#'
#' @format Include the heritability, the model used (models with squared brackets use genotypic information), the number of varieties, the number of observations, the Kendall metric and the top20 metric
#'
#' @keywords datasets
#'
#' @examples
#' data(resOutlier)
#' summary(resOutlier)

"resOutlier"
