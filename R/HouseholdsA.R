#' @name HouseholdsA
#' @aliases HouseholdsA
#' @docType data
#' @title Database of household expenses for frame A
#' 
#' @description This dataset contains some variables regarding household expenses for a sample of 105 households selected from a list of landline phones (let say, frame A) in a particular city in a specific month.
#' @usage HouseholdsA
#' @details The sample, of size \eqn{n_A = 105}, has been drawn from a population of \eqn{N_A = 1735} households with landline phone. \eqn{N_ab = 601} of these households have, also, mobile phone.
#'  On the other hand, frame totals for auxiliary variables in this frame are \eqn{X_{Income}^A = 4300260} and \eqn{X_{Taxes}^A = 215577}.
#' @format
#' \describe{
#'    \item{Domain}{A string indicating the domain each household belongs to. Possible values are "a" if household belongs to domain a or "ab" if household belongs to overlap domain.}
#'    \item{Feeding}{Feeding expenses (in euros) at the househould}
#'    \item{Clothing}{Clothing expenses (in euros) at the household}
#'    \item{Leisure}{Leisure expenses (in euros) at the household}
#'    \item{Income}{Household income (in euros). Values for this variable are only available for households included in frame A.}
#'    \item{Taxes}{Household municipal taxes (in euros) paid. Values for this variable are only available for households included in frame A.}
#'    \item{Metres2}{Square meters of the house. Values for this variable are only available for households included in frame B.}
#'    \item{Size}{Household size. Values for this variable are only available for households included in frame B.}
#'    \item{ProbA}{First order inclusion probability in frame A.}
#'    \item{ProbB}{First order inclusion probability in frame B.}
#' }
#' @examples
#' data(HouseholdsA)
#' attach(HouseholdsA)
#' #Let perform a brief descriptive analysis for the three main variables
#' param <- data.frame(Feeding, Clothing, Leisure)
#' summary (param)
#' hist (Feeding)
#' hist (Clothing)
#' hist (Leisure)
#' @keywords datasets
NULL