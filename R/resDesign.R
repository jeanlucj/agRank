#' Simulation results for comparison of experimental designs
#'
#' In this simulation,
#' we only include LM and BT with genotypic information.
#' We run the procedure 100 times with h^2 = 0.2, 0.5, 0.8 with m = 240 and n = 180, 240, 300, 360.
#' 17 Experimental designs are included, namely AD, RD, GD1~GD9, KM1~KM6. Please see the vignette for detailed
#' description of each experimental design.
#'
#' @docType data
#'
#' @usage data(resDesign)
#'
#' @format Include the heritability, the design used, the model used (models with squared brackets use genotypic information), the number of varieties, the number of observations, the Kendall metric and the top20 metric
#'
#' @keywords datasets
#'
#' @examples
#' data(resDesign)
#' summary(resDesign)

"resDesign"
