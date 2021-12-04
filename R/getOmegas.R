#' Return omegas 
#'
#' @param k category index
#' @param omega list of omegas
#'
#' @return vector of omegas for category k 
#'

getOmegas <- function(k, omega){
  
  n = length(omega)
  
  unlist(sapply(1:n, FUN = function(i){
    if(!is.na(omega[[i]][k])) return(omega[[i]][k])
    else return(NULL)
  }))
}