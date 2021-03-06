% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Hartley.R
\name{Hartley}
\alias{Hartley}
\title{Hartley estimator}
\usage{
Hartley(ysA, ysB, pi_A, pi_B, domains_A, domains_B, conf_level = NULL)
}
\arguments{
\item{ysA}{A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable of interest from \eqn{s_A}.}

\item{ysB}{A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable of interest from \eqn{s_B}.}

\item{pi_A}{A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.}

\item{pi_B}{A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.}

\item{domains_A}{A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".}

\item{domains_B}{A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".}

\item{conf_level}{(Optional) A numeric value indicating the confidence level for the confidence intervals.}
}
\value{
\code{Hartley} returns an object of class "EstimatorDF" which is a list with, at least, the following components:
 \item{Call}{the matched call.}
 \item{Est}{total and mean estimation for main variable(s).}
 \item{VarEst}{variance estimation for main variable(s).}
 If parameter \code{conf_level} is different from \code{NULL}, object includes component
 \item{ConfInt}{total and mean estimation and confidence intervals for main variables(s).}
 In addition, components \code{TotDomEst} and \code{MeanDomEst} are available when estimator is based on estimators of the domains. Component \code{Param} shows value of parameters involded in calculation of the estimator (if any).
 By default, only \code{Est} component (or \code{ConfInt} component, if parameter \code{conf_level} is different from \code{NULL}) is shown. It is possible to access to all the components of the objects by using function \code{summary}.
}
\description{
Produces estimates for population totals and means using Hartley estimator from survey data obtained
 from a dual frame sampling design. Confidence intervals are also computed, if required.
}
\details{
Hartley estimator of population total is given by
 \deqn{\hat{Y}_H = \hat{Y}_a^A + \hat{\theta}\hat{Y}_{ab}^A + (1 - \hat{\theta})\hat{Y}_{ab}^B + \hat{Y}_b^B}
 where \eqn{\hat{\theta} \in [0, 1]}. Optimum value for \eqn{\hat{\theta}} to minimize variance of the estimator is
 \deqn{\hat{\theta}_{opt} = \frac{\hat{V}(\hat{Y}_{ab}^B) + \widehat{Cov}(\hat{Y}_b^B, \hat{Y}_{ab}^B) - \widehat{Cov}(\hat{Y}_a^A, \hat{Y}_{ab}^A)}{\hat{V}(\hat{Y}_{ab}^A) + \hat{V}(\hat{Y}_{ab}^B)}}
 Taking into account the independence between \eqn{s_A} and \eqn{s_B}, an estimator for the variance of the Hartley estimator can be obtained as follows:
 \deqn{\hat{V}(\hat{Y}_H) = \hat{V}(\hat{Y}_a^A + \hat{\theta}\hat{Y}_{ab}^A) + \hat{V}((1 - \hat{\theta})\hat{Y}_{ab}^B + \hat{Y}_b^B)}
 If both first and second order probabilities are known, variances and covariances involved in calculation of \eqn{\hat{\theta}_{opt}} and \eqn{\hat{V}(\hat{Y}_H)} are estimated using functions \code{VarHT} and \code{CovHT}, respectively. If
 only first order probabilities are known, variances are estimated using Deville's method and covariances are estimated using following expression
 \deqn{\widehat{Cov}(\hat{X}, \hat{Y}) = \frac{\hat{V}(X + Y) - \hat{V}(X) - \hat{V}(Y)}{2}}
}
\examples{
data(DatA)
data(DatB)
data(PiklA)
data(PiklB)

#Let calculate Hartley estimator for variable Feeding
Hartley(DatA$Feed, DatB$Feed, PiklA, PiklB, DatA$Domain, DatB$Domain)

#Now, let calculate Hartley estimator and a 90\% confidence interval
#for variable Leisure, considering only first order inclusion probabilities
Hartley(DatA$Lei, DatB$Lei, DatA$ProbA, DatB$ProbB, DatA$Domain, 
DatB$Domain, 0.90)
}
\references{
Hartley, H. O. (1962)
 \emph{Multiple Frames Surveys.}
 Proceedings of the American Statistical Association, Social Statistics Sections, 203 - 206.

Hartley, H. O. (1974) 
 \emph{Multiple frame methodology and selected applications.}
 Sankhya C, Vol. 36, 99 - 118.
}
\seealso{
\code{\link{JackHartley}}
}

