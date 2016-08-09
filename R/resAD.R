#' Simulation results for models with adherence parameters (except LM)
#'
#'We do a simulation to compare the models.
#'Namely we include LM, TH, BT and PL with and without genotypic information.
#'All models except LM can accommodate the adherence parameters.
#'We break this simulation into two parts: the pilot simulation and the main simulation.
#'In the pilot simulation, we run the simulation 100 times with h^2 = 0.2, 0.5, 0.8,
#'for nvar = 15, 30, 60 and nobs = 30, 60, 120.
#'In the main simulation, we run the simulation 100 times with
#'h^2 = 0.2, 0.5, 0.8, for nvar = 120, 240, 360 and nobs = 180, 240, 300, 360.
#'We use the alpha design to allocate varieties to farmers.
#'
#' @docType data
#'
#' @usage data(resAD)
#'
#' @format Include the heritability, the component (main or pilot simulation), the model used (models with squared brackets use genotypic information), the number of varieties, the number of observations, the Kendall metric and the top20 metric
#'
#' @keywords datasets
#'
#' @examples
#' data(resAD)
#' summary(resAD)

"resAD"
