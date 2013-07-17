#' @name Hartley
#' @aliases Hartley
#' @title Hartley estimator
#' 
#' @description Produces estimates for population totals and means using Hartley estimator from survey data obtained
#'  from a dual frame sampling design. Confidence intervals are also computed, if required.
#' 
#' @usage Hartley(ysA, ysB, pi_A, pi_B, domains_A, domains_B, conf_level = NULL)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals.
#' @details Hartley estimator of population total is given by
#'  \deqn{\hat{Y}_H = \hat{Y}_a^A + \hat{\theta}\hat{Y}_{ab}^A + (1 - \hat{\theta})\hat{Y}_{ab}^B + \hat{Y}_b^B}
#'  where \eqn{\hat{\theta} \in [0, 1]}. Optimum value for \eqn{\hat{\theta}} to minimize variance of the estimator is
#'  \deqn{\hat{\theta}_{opt} = \frac{\hat{V}(\hat{Y}_{ab}^B) + \widehat{Cov}(\hat{Y}_b^B, \hat{Y}_{ab}^B) - \widehat{Cov}(\hat{Y}_a^A, \hat{Y}_{ab}^A)}{\hat{V}(\hat{Y}_{ab}^A) + \hat{V}(\hat{Y}_{ab}^B)}}
#'  Taking into account the independence between \eqn{s_A} and \eqn{s_B}, an estimator for the variance of the Hartley estimator can be obtained as follows:
#'  \deqn{\hat{V}(\hat{Y}_H) = \hat{V}(\hat{Y}_a^A + \hat{\theta}\hat{Y}_{ab}^A) + \hat{V}((1 - \hat{\theta})\hat{Y}_{ab}^B + \hat{Y}_b^B)}
#'  If both first and second order probabilities are known, variances and covariances involved in calculation of \eqn{\hat{\theta}_{opt}} and \eqn{\hat{V}(\hat{Y}_H)} are estimated using functions \code{VarHT} and \code{CovHT}, respectively. If
#'  only first order probabilities are known, variances are estimated using Deville's method and covariances are estimated using following expression
#'  \deqn{\widehat{Cov}(\hat{X}, \hat{Y}) = \frac{\hat{V}(X + Y) - \hat{V}(X) - \hat{V}(Y)}{2}} 
#' @return A numeric matrix containing estimations of population total and population mean for considered values.
#' @references Hartley, H. O. (1962)
#'  \emph{Multiple Frames Surveys.}
#'  Proceedings of the American Statistical Association, Social Statistics Sections, 203 - 206.
#' @references Hartley, H. O. (1974) 
#'  \emph{Multiple frame methodology and selected applications.}
#'  Sankhya C, Vol. 36, 99 - 118.
#' @seealso \code{\link{JackHartley}}
#' @examples
#' data(HouseholdsA)
#' dataA <- attach(HouseholdsA)
#' detach (HouseholdsA)
#' data(HouseholdsB)
#' dataB <- attach(HouseholdsB)
#' detach(HouseholdsB)
#' data(PiklA)
#' data(PiklB)
#' 
#' #Let calculate Hartley estimator for variable Feeding
#' Hartley(dataA$Feeding, dataB$Feeding, PiklA, PiklB, dataA$Domain, dataB$Domain)
#' 
#' #Now, let calculate Hartley estimator and a 90% confidence interval
#' #for variable Leisure, considering only first order inclusion probabilities
#' Hartley(dataA$Leisure, dataB$Leisure, dataA$ProbA, dataB$ProbB, dataA$Domain, 
#' dataB$Domain, 0.90)
#' @export
Hartley = function (ysA, ysB, pi_A, pi_B, domains_A, domains_B, conf_level = NULL)
{	
	cnames <- names(ysA)
	ysA <- as.matrix(ysA)
	ysB <- as.matrix(ysB)
	pi_A <- as.matrix(pi_A)
	pi_B <- as.matrix(pi_B)
	
	if (any(is.na(ysA)))
		stop("There are missing values in sample from frame A.")
	if (any(is.na(ysB)))
		stop("There are missing values in sample from frame B.")
	if (any(is.na(pi_A)))
		stop("There are missing values in inclusion probabilities from frame A.")
	if (any(is.na(pi_B)))
		stop("There are missing values in inclusion probabilities from frame B.")
	if (any(is.na(domains_A)))
		stop("There are missing values in domains from frame A.")
	if (any(is.na(domains_B)))
		stop("There are missing values in domains from frame B.")
	if (nrow(ysA) != nrow(pi_A) | nrow(ysA) != length(domains_A) | length(domains_A) != nrow(pi_A))
		stop("Arguments from frame A have different sizes.")
	if (nrow(ysB) != nrow(pi_B) | nrow(ysB) != length(domains_B) | length(domains_B) != nrow(pi_B))
		stop("Arguments from frame B have different sizes.")
	if (ncol(ysA) != ncol(ysB))
		stop("Number of variables does not match.")
	if (length(which(domains_A == "a")) + length(which(domains_A == "ab")) != length(domains_A))
		stop("Domains from frame A are not correct.")
	if (length(which(domains_B == "b")) + length(which(domains_B == "ba")) != length(domains_B))
		stop("Domains from frame B are not correct.")

	if (is.null(conf_level)) {
		r <- 2
		rnames <- c("Total", "Mean")
	}
	else {
		r <- 6
		rnames <- c("Total", "Upper End", "Lower End", "Mean", "Upper End", "Lower End")
	}

	n_A <- nrow(ysA)
	n_B <- nrow(ysB)
	c <- ncol(ysA)
	results <- matrix(NA, nrow = r, ncol = c)
	rownames(results) <- rnames
	colnames(results) <- cnames

	ysA <- cbind(rep(1, n_A), ysA)
	ysB <- cbind(rep(1, n_B), ysB)
	c1 <- ncol(ysA)
	
	if (!is.null(dim(drop(pi_A))) & !is.null(dim(drop(pi_B)))) {

		if (nrow(pi_A) != ncol(pi_A))
			stop("Matrix of inclusion probabilities from frame A is not square.")
		if (nrow(pi_B) != ncol(pi_B))
			stop("Matrix of inclusion probabilities from frame B is not square.")

		for (k in 1:c1) {

			data_a_A <- Domains(ysA[,k], domains_A, "a")
			data_ab_A <- Domains(ysA[,k], domains_A, "ab")
			data_b_B <- Domains(ysB[,k], domains_B, "b")
			data_ab_B <- Domains(ysB[,k], domains_B, "ba")

			Yhat_a_A <- HT (data_a_A, diag(pi_A))
			Yhat_ab_A <- HT (data_ab_A, diag(pi_A))
			Yhat_b_B <- HT (data_b_B, diag(pi_B))
			Yhat_ab_B <- HT (data_ab_B, diag(pi_B))

			Vhat_Yhat_a_A <- VarHT (data_a_A, pi_A)
			Vhat_Yhat_ab_A <- VarHT (data_ab_A, pi_A)
			Vhat_Yhat_b_B <- VarHT (data_b_B, pi_B)
			Vhat_Yhat_ab_B <- VarHT (data_ab_B, pi_B)
			Covhat_Yhat_a_A_Yhat_ab_A <- CovHT (data_a_A, data_ab_A, pi_A)
			Covhat_Yhat_b_B_Yhat_ab_B <- CovHT (data_b_B, data_ab_B, pi_B)
			    
			theta_opt <- (Vhat_Yhat_ab_B + Covhat_Yhat_b_B_Yhat_ab_B - Covhat_Yhat_a_A_Yhat_ab_A) / (Vhat_Yhat_ab_A + Vhat_Yhat_ab_B)
			if (theta_opt > 1 | theta_opt < 0) {
				warning ("Optimal value of theta is smaller than 0 or greater than 1. Estimation V(Y_ab^B) / (V(Y_ab^A) + V(Y_ab^B)) has been used instead.")
				theta_opt <- Vhat_Yhat_ab_B / (Vhat_Yhat_ab_A + Vhat_Yhat_ab_B)
			}

			if (k == 1)
				size_estimation <- Yhat_a_A + theta_opt * Yhat_ab_A + (1 - theta_opt) * Yhat_ab_B + Yhat_b_B
			else
				total_estimation <- Yhat_a_A + theta_opt * Yhat_ab_A + (1 - theta_opt) * Yhat_ab_B + Yhat_b_B

			if (k > 1) {
				mean_estimation <- total_estimation / size_estimation
				results[,k-1] <- c(total_estimation, mean_estimation)

				if (!is.null(conf_level)) {

					Vhat_Yhat_Hartley <- Vhat_Yhat_a_A + theta_opt^2 * Vhat_Yhat_ab_A + (1 - theta_opt)^2 * Vhat_Yhat_ab_B + Vhat_Yhat_b_B + 2 * theta_opt * Covhat_Yhat_a_A_Yhat_ab_A + 2 * (1 - theta_opt) * Covhat_Yhat_b_B_Yhat_ab_B
      					total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_Hartley)
					total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_Hartley)
					mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * Vhat_Yhat_Hartley)
					mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * Vhat_Yhat_Hartley)
					results[,k-1] <- c(total_estimation, total_upper, total_lower, mean_estimation, mean_upper, mean_lower)
				}
			}
		}
	}
	else {
		
		if (is.null(dim(drop(pi_A))) & is.null(dim(drop(pi_B)))){

			for (k in 1:c1) {

				data_a_A <- Domains(ysA[,k], domains_A, "a")
				data_ab_A <- Domains(ysA[,k], domains_A, "ab")
				data_b_B <- Domains(ysB[,k], domains_B, "b")
				data_ab_B <- Domains(ysB[,k], domains_B, "ba")

				Yhat_a_A <- HT (data_a_A, pi_A)
				Yhat_ab_A <- HT (data_ab_A, pi_A)
				Yhat_b_B <- HT (data_b_B, pi_B)
				Yhat_ab_B <- HT (data_ab_B, pi_B)

				Vhat_Yhat_a_A <- varest (data_a_A, pik = pi_A)
				Vhat_Yhat_ab_A <- varest (data_ab_A, pik = pi_A)
				Vhat_Yhat_b_B <- varest (data_b_B, pik = pi_B)
				Vhat_Yhat_ab_B <- varest (data_ab_B, pik = pi_B)
				Covhat_Yhat_a_A_Yhat_ab_A <- (varest (Ys = ysA[,k], pik = pi_A) - Vhat_Yhat_a_A - Vhat_Yhat_ab_A)/2
				Covhat_Yhat_b_B_Yhat_ab_B <- (varest (Ys = ysB[,k], pik = pi_B) - Vhat_Yhat_b_B - Vhat_Yhat_ab_B)/2

				theta_opt <- (Vhat_Yhat_ab_B + Covhat_Yhat_b_B_Yhat_ab_B - Covhat_Yhat_a_A_Yhat_ab_A) / (Vhat_Yhat_ab_A + Vhat_Yhat_ab_B)
				if (theta_opt > 1 | theta_opt < 0) {
					warning ("Optimal value of theta is smaller than 0 or greater than 1. Estimation V(Y_ab^B) / (V(Y_ab^A) + V(Y_ab^B)) has been used instead.")
					theta_opt <- Vhat_Yhat_ab_B / (Vhat_Yhat_ab_A + Vhat_Yhat_ab_B)
				}

				if (k == 1)
					size_estimation <- Yhat_a_A + theta_opt * Yhat_ab_A + (1 - theta_opt) * Yhat_ab_B + Yhat_b_B
				else
					total_estimation <- Yhat_a_A + theta_opt * Yhat_ab_A + (1 - theta_opt) * Yhat_ab_B + Yhat_b_B

				if (k > 1) {
					mean_estimation <- total_estimation/size_estimation
					results[,k-1] <- c(total_estimation, mean_estimation)

					if (!is.null(conf_level)) {
			
						Vhat_Yhat_Hartley <- Vhat_Yhat_a_A + theta_opt^2 * Vhat_Yhat_ab_A + (1 - theta_opt)^2 * Vhat_Yhat_ab_B + Vhat_Yhat_b_B + 2 * theta_opt * Covhat_Yhat_a_A_Yhat_ab_A + 2 * (1 - theta_opt) * Covhat_Yhat_b_B_Yhat_ab_B
						total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_Hartley)
						total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_Hartley)
						mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * Vhat_Yhat_Hartley)
						mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * Vhat_Yhat_Hartley)
						results[,k-1] <- c(total_estimation, total_upper, total_lower, mean_estimation, mean_upper, mean_lower)
					}
				}	
			}
		}
		else
			stop("Invalid option: Probability vector in one frame and probability matrix in the other frame. Type of probabilities structures must match.")
	}
	return (results)
}