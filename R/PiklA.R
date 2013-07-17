#' @name PiklA
#' @aliases PiklA
#' @docType data
#' @title Matrix of inclusion probabilities for units selected in sample from frame A
#' 
#' @description This dataset consists of a square matrix of dimension 105 with the first and second order inclusion probabilities
#'  for the units included in \eqn{s_A}.
#' @usage PiklA
#' 
#' @examples 
#' data(PiklA)
#' #Let choose the submatrix of inclusion probabilities for the first 5 units sA.
#' PiklA[1:5, 1:5]
#' #Now, let select only the first order inclusion probabilities
#' diag(PiklA)
#' 
#' @keywords datasets
NULL