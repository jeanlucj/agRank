#' Marker matrix used in the simulations
#'
#' The genotype data are wheat marker sequences from Kansas State Genotyping Lab (KSG),
#' available in the Triticeae Toolbox.
#' Marker data are filtered by maximum missing data = 8\% and minimum MAF = 5\%,
#' and we get 368 lines with 9928 markers for each line. The marker matrix is coded in (-1, 0, 1).
#'
#' @docType data
#'
#' @usage data(marker)
#'
#'
#' @keywords datasets
#'
#' @source \url{https://triticeaetoolbox.org/wheat/display_genotype.php?function=typeData&mm=8&mmaf=5&trial_code=2014_TCAP_ABB_SRW_Mid}
#'
#' @examples
#' data(marker)
#' dim(marker) #368 lines, 9928 loci

"marker"
