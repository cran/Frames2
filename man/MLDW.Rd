% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MLDW.R
\name{MLDW}
\alias{MLDW}
\title{Multinomial logistic estimator under dual frame approach with auxiliary information from the whole population}
\usage{
MLDW (ysA, ysB, pik_A, pik_B, domains_A, domains_B, xsA, xsB, x, ind_sam, 
 conf_level = NULL)
}
\arguments{
\item{ysA}{A data frame containing information about one or more factors, each one of dimension \eqn{n_A}, collected from \eqn{s_A}.}

\item{ysB}{A data frame containing information about one or more factors, each one of dimension \eqn{n_B}, collected from \eqn{s_B}.}

\item{pik_A}{A numeric vector of length \eqn{n_A} containing first order inclusion probabilities for units included in \eqn{s_A}.}

\item{pik_B}{A numeric vector of length \eqn{n_B} containing first order inclusion probabilities for units included in \eqn{s_B}.}

\item{domains_A}{A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".}

\item{domains_B}{A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".}

\item{xsA}{A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information in frame A for units included in \eqn{s_A}.}

\item{xsB}{A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information in frame B for units included in \eqn{s_B}.}

\item{x}{A numeric vector or length \eqn{N} or a numeric matrix or data frame of dimensions \eqn{N} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information for every unit in the population.}

\item{ind_sam}{A numeric vector of length \eqn{n = n_A + n_B} containing the identificators of units of the population (from 1 to \eqn{N}) that belongs to \eqn{s_A} or \eqn{s_B}}

\item{conf_level}{(Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.}
}
\value{
\code{MLDW} returns an object of class "MultEstimatorDF" which is a list with, at least, the following components:
 \item{Call}{the matched call.}
 \item{Est}{class frequencies and proportions estimations for main variable(s).}
}
\description{
Produces estimates for class totals and proportions using multinomial logistic regression from survey data obtained
 from a dual frame sampling design using a dual frame model assisted approach. Confidence intervals are also computed, if required.
}
\details{
Multinomial logistic estimator in dual frame using auxiliary information from the whole population for a proportion is given by
 \deqn{\hat{P}_{MLi}^{DW} = \frac{1}{N} (\sum_{k \in U} p_{ki}^{\circ} + \sum_{k \in s} {d}_k^{\circ} (z_{ki} - p_{ki}^{\circ})) \hspace{0.3cm} i = 1,...,m}
 with \eqn{m} the number of categories of the response variable, \eqn{z_i} the indicator variable for the i-th category of the response variable,
 \eqn{d_k^{\circ} =\left\{\begin{array}{lcc}
 d_k^A & \textrm{if } k \in a\\
 \eta d_k^A & \textrm{if } k \in ab\\
 (1 - \eta) d_k^B & \textrm{if } k \in ba \\
 d_k^B & \textrm{if } k \in b
 \end{array}
 \right.}
 with \eqn{\eta \in (0,1)} and \deqn{p_{ki}^\circ = \frac{exp(x_k^{'}\beta_i^{\circ})}{\sum_{r=1}^m exp(x_k^{'}\beta_r^{\circ})},}
being \eqn{\beta_i^{\circ}} the maximum likelihood parameters of the multinomial logistic model considering the weights \eqn{d^{\circ}}.
}
\examples{
data(DatMA)
data(DatMB)
data(DatPopM) 

IndSample <- c(DatMA$Id_Pop, DatMB$Id_Pop)
#Let calculate proportions of categories of variable Prog using MLDW estimator
#using Read as auxiliary variable
MLDW(DatMA$Prog, DatMB$Prog, DatMA$ProbA, DatMB$ProbB, DatMA$Domain, DatMB$Domain, 
DatMA$Read, DatMB$Read, DatPopM$Read, IndSample)

#Let obtain 95\% confidence intervals together with the estimations
MLDW(DatMA$Prog, DatMB$Prog, DatMA$ProbA, DatMB$ProbB, DatMA$Domain, DatMB$Domain, 
DatMA$Read, DatMB$Read, DatPopM$Read, IndSample, 0.95)
}
\references{
Molina, D., Rueda, M., Arcos, A. and Ranalli, M. G. (2015)
 \emph{Multinomial logistic estimation in dual frame surveys}
 Statistics and Operations Research Transactions (SORT). To be printed.

Lehtonen, R. and Veijanen, A. (1998) 
 \emph{On multinomial logistic generalizaed regression estimators}
 Technical report 22, Department of Statistics, University of Jyvaskyla.
}
\seealso{
\code{\link{JackMLDW}}
}

