#' Allocation of varieties to farmers
#'
#' Generate experimental designs to allocate vartieties to farmers.
#' two experimental designs are included, "alpha" and "random".
#'
#' @param nVarieties the number of varieties
#' @param nFarmers the number of farmers (observers), must be greater than or equal to nVarieties / 3
#' @param method must be one of "alpha" or "random"
#'
#' @return an nFarmers * 3 matrix, with each row being the label of varieties assigned to that farmer
#'
#' @examples
#' #use alpha design
#' expDesign(120, 245, method = 'alpha')
#'
#' @export
#'
expDesign <- function(nVarieties, nFarmers, method=c("alpha", "random")){

  switch(method,
         alpha = {
           #alpha design
           library(agricolae)
           nRep <- nFarmers*3 / nVarieties
           alphAssign <- NULL
           for (i in 1:ceiling(nRep / 3)){
             sketch <- design.alpha(1:nVarieties, 3, 3)$sketch
             alphAssign <- rbind(alphAssign, sketch[[1]], sketch[[2]], sketch[[3]])
           }
           alphAssign <- alphAssign[1:nFarmers,]
           alphAssign <- matrix(as.numeric(alphAssign), nFarmers)
         },#END alpha method

         random = {
           nPlots <- 3 * nFarmers
           randAssign <- matrix(sample(rep(1:nVarieties, length.out=nPlots)), nrow=nFarmers)
           # Make sure randAssign does not have the same variety twice for one farmer
           nUni <- apply(randAssign, 1, function(vec) length(unique(vec))) # if 3 then good
           for (farmer in which(nUni < 3)){
             # find other farmers that don't have the problem
             while (length(unique(randAssign[farmer,])) < 3){
               nOfVar <- table(randAssign[farmer,])
               # Which variety is doubled?
               problemVar <- as.numeric(names(nOfVar)[nOfVar > 1])
               # Which farmers were not assigned that variety?
               posSwap <- apply(randAssign, 1, function(vec) sum(vec == problemVar) == 0)
               # Swap with one of those farmers
               swapFarmer <- sample(which(posSwap), 1)
               randAssign[farmer, which(randAssign[farmer,] == problemVar)[1]] <- randAssign[swapFarmer, 1]
               randAssign[swapFarmer, 1] <- problemVar
             }
           }
           randAssign
         }#END random method
  )#END switch
}
