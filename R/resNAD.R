#' Simulation results for models without adherence parameters
#'
#'We do a simulation to see the performance of models
#'when the adherence parameters are not included.
#'Models without adherence parameters are simply recovered by taking adherences for all farmers to be 1.
#'We include LM, TH, BT and PL with and without genetypic information in this simulation.
#'We run the procedure 100 times with h^2 = 0.2, 0.5, 0.8 with nvar = 240 and nobs = 180, 240, 300, 360.
#'We use the alpha design to allocate varieties to farmers.
#'
#'
#' @docType data
#'
#'
#' @usage data(resNAD)
#'
#' @format Include the heritability, the model used (models with squared brackets use genotypic information), the number of varieties, the number of observations, the Kendall metric and the top20 metric
#'
#' @keywords datasets
#'
#' @examples
#' data(resNAD)
#' summary(resNAD)

"resNAD"
