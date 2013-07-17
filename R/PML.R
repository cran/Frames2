#' @name PML
#' @aliases PML
#' @title Pseudo Maximum Likelihood estimator
#' 
#' @description Produces estimates for population totals and means using PML estimator from survey data obtained
#'  from a dual frame sampling design. Confidence intervals are also computed, if required.
#' 
#' @usage PML(ysA, ysB, pi_A, pi_B, domains_A, domains_B, N_A, N_B, conf_level = NULL)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param N_A A numeric value indicating the size of frame A
#' @param N_B A numeric value indicating the size of frame B
#' @param conf_level (Optional) A numeric value indicating the confidence level for the confidence intervals, if desired.
#' @details Pseudo Maximum Likelihood estimator of population total is given by
#'  \deqn{\hat{Y}_{PML}(\hat{\theta}) = \frac{N_A - \hat{N}_{ab,PML}}{\hat{N}_a}\hat{Y}_a^A + \frac{N_B - \hat{N}_{ab,PML}}{\hat{N}_b}\hat{Y}_b^B + \frac{\hat{N}_{ab,PML}}{\hat{\theta}\hat{N}_{ab}^A + (1 - \hat{\theta})\hat{N}_{ab}^B}[\hat{\theta}\hat{Y}_{ab}^A + (1 - \hat{\theta})\hat{Y}_{ab}^B]}
#'  where \eqn{\hat{\theta} \in [0, 1]} and \eqn{\hat{N}_{ab,PML}} is the smaller of the roots of the quadratic equation
#'  \deqn{[\hat{\theta}/N_B + (1 - \hat{\theta})/N_A]x^2 - [1 + \hat{\theta}\hat{N}_{ab}^A/N_B + (1 - \hat{\theta})\hat{N}_{ab}^B/N_A]x + \hat{\theta}\hat{N}_{ab}^A + (1 - \hat{\theta})\hat{N}_{ab}^B=0.} Optimal value for \eqn{\hat{\theta}} is \eqn{\frac{\hat{N}_aN_B\hat{V}(\hat{N}_{ab}^B)}{\hat{N}_aN_B\hat{V}(\hat{N}_{ab}^B) + \hat{N}_bN_A\hat{V}(\hat{N}_{ab}^A)}} 
#' @return A numeric matrix containing estimations of population total and population mean for considered values.
#' @references Skinner, C. J. and Rao, J. N. K. (1996)
#'  \emph{Estimation in Dual Frame Surveys with Complex Designs}. Journal of the American Statistical Association, Vol. 91, 433, 349 - 356.
#' @examples
#' data(HouseholdsA)
#' dataA <- attach(HouseholdsA)
#' detach(HouseholdsA)
#' data(HouseholdsB)
#' dataB <- attach(HouseholdsB)
#' detach(HouseholdsB)
#' data(PiklA)
#' data(PiklB)
#'
#' #Let calculate Pseudo Maximum Likelihood estimator for population total for variable Clothing
#' PML(dataA$Clothing, dataB$Clothing, PiklA, PiklB, dataA$Domain, dataB$Domain, 
#' N_A = 1735, N_B = 1191)
#' 
#' #Now, let calculate Pseudo Maximum Likelihood estimator for population total for variable
#' #Feeding, using first order inclusion probabilities
#'
#' PML(dataA$Feeding, dataB$Feeding, dataA$ProbA, dataB$ProbB, dataA$Domain, dataB$Domain, 
#' N_A = 1735, N_B = 1191)
#'
#' #Finally, let calculate Pseudo Maximum Likelihood estimator and a 90% confidence interval for 
#' #population total for variable Leisure
#' PML(dataA$Leisure, dataB$Leisure, PiklA, PiklB, dataA$Domain, dataB$Domain, 
#' N_A = 1735, N_B = 1191, 0.90)
#' @export
PML = function (ysA, ysB, pi_A, pi_B, domains_A, domains_B, N_A, N_B, conf_level = NULL)
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

	ones_a_A <- Domains (rep (1, n_A), domains_A, "a")
	ones_ab_A <- Domains (rep (1, n_A), domains_A, "ab")
	ones_b_B <- Domains (rep (1, n_B), domains_B, "b")
	ones_ab_B <- Domains (rep (1, n_B), domains_B, "ba")
	
	if (!is.null(dim(drop(pi_A))) & !is.null(dim(drop(pi_B)))) {

		if (nrow(pi_A) != ncol(pi_A))
			stop("Matrix of inclusion probabilities from frame A is not square.")
		if (nrow(pi_B) != ncol(pi_B))
			stop("Matrix of inclusion probabilities from frame B is not square.")
	
		Nhat_a_A <- HT (ones_a_A, diag(pi_A))
		Nhat_b_B <- HT (ones_b_B, diag(pi_B))
		Nhat_ab_A <- HT (ones_ab_A, diag(pi_A))
		Nhat_ab_B <- HT (ones_ab_B, diag(pi_B))
		Vhat_Nhat_ab_B <- VarHT (ones_ab_B, pi_B)
		Vhat_Nhat_ab_A <- VarHT (ones_ab_A, pi_A)

		for (k in 1:c1) {
	
			data_a_A <- Domains (ysA[,k], domains_A, "a")
			data_ab_A <- Domains (ysA[,k], domains_A, "ab")
			data_b_B <- Domains (ysB[,k], domains_B, "b")
			data_ab_B <- Domains (ysB[,k], domains_B, "ba")

			Yhat_a_A <- HT(data_a_A, diag(pi_A))
			Yhat_ab_A <- HT (data_ab_A, diag(pi_A))
			Yhat_b_B <- HT (data_b_B, diag(pi_B))
			Yhat_ab_B <- HT (data_ab_B, diag(pi_B))

			theta_opt <- (Nhat_a_A * N_B * Vhat_Nhat_ab_B) / (Nhat_a_A * N_B * Vhat_Nhat_ab_B + Nhat_b_B * N_A * Vhat_Nhat_ab_A)
			if (is.nan(theta_opt) | theta_opt == 0 | theta_opt == 1){ #Alg?n dominio es vac?o (o los dos). Entonces calculamos los efectos del dise?o para obtener los tama?os muestrales efectivos
		
				SRSPikls_A <- matrix((n_A * (n_A - 1)) / (N_A * (N_A - 1)), n_A, n_A)
				diag(SRSPikls_A) <- rep(n_A / N_A, n_A)
				SRSPikls_B <- matrix((n_B * (n_B - 1)) / (N_B * (N_B - 1)), n_B, n_B)
				diag(SRSPikls_B) <- rep(n_B / N_B, n_B)
	
				d_A <- Vhat_Nhat_ab_A / VarHT (ones_ab_A, SRSPikls_A) 
				d_B <- Vhat_Nhat_ab_B / VarHT (ones_ab_B, SRSPikls_B)

				eff_n_A <- n_A / d_A
				eff_n_B <- n_B / d_B
		
				theta_opt <- eff_n_A * N_B / (eff_n_A * N_B + eff_n_B * N_A)
			}
			term_a <- theta_opt / N_B + (1 - theta_opt) / N_A
			term_b <- -(1 + theta_opt * Nhat_ab_A / N_B + (1 - theta_opt) * Nhat_ab_B / N_A)
			term_c <- theta_opt * Nhat_ab_A + (1 - theta_opt) * Nhat_ab_B
			Nhat_ab_PML_theta_opt <- (- term_b - sqrt(term_b * term_b - 4 * term_a * term_c)) / (2 * term_a)

			muhat_ab <- (n_A / N_A * Yhat_ab_A + n_B / N_B * Yhat_ab_B) / (n_A / N_A * Nhat_ab_A + n_B / N_B * Nhat_ab_B)
			muhat_a <- Yhat_a_A / Nhat_a_A
			muhat_b <- Yhat_b_B / Nhat_b_B

			if (k == 1)
				size_estimation <- sum((N_A - Nhat_ab_PML_theta_opt) * muhat_a, (N_B - Nhat_ab_PML_theta_opt) * muhat_b, Nhat_ab_PML_theta_opt * muhat_ab, rm.na = TRUE)
			else
				total_estimation <- sum((N_A - Nhat_ab_PML_theta_opt) * muhat_a, (N_B - Nhat_ab_PML_theta_opt) * muhat_b, Nhat_ab_PML_theta_opt * muhat_ab, rm.na = TRUE) 
		
			if (k > 1) {
				mean_estimation <- total_estimation / size_estimation
				results[,k-1] <- c(total_estimation, mean_estimation)

				if (!is.null(conf_level)) {
			
					lambdahat <- sum(muhat_ab, - muhat_a, - muhat_b, na.rm = TRUE)
					phihat <- n_A * Nhat_b_B / (n_A * Nhat_b_B + n_B * Nhat_a_A)
					if (is.nan(phihat))
						phihat <- 0

					zA1 <- cbind((ysA[,k] - muhat_a) * ones_a_A, (theta_opt * (ysA[,k] - muhat_ab) + lambdahat * phihat) * ones_ab_A)
					zB1 <- cbind((ysB[,k] - muhat_b) * ones_b_B, ((1 - theta_opt) * (ysB[,k] - muhat_ab) + lambdahat * (1 - phihat)) * ones_ab_B)
					zA <- rowSums(zA1, na.rm = TRUE)
					zB <- rowSums(zB1, na.rm = TRUE)

					Vhat_Yhat_PML <- VarHT (zA, pi_A) + VarHT (zB, pi_B)
					total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_PML)
					total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_PML)
					mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * Vhat_Yhat_PML)
					mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * Vhat_Yhat_PML)
					results[,k-1] <- c(total_estimation, total_upper, total_lower, mean_estimation, mean_upper, mean_lower)
				
				}
			}
		}
	}
	else {

		if(is.null(dim(drop(pi_A))) & is.null(dim(drop(pi_B)))){
			
			Nhat_a_A <- HT (ones_a_A, pi_A)
			Nhat_b_B <- HT (ones_b_B, pi_B)
			Nhat_ab_A <- HT (ones_ab_A, pi_A)
			Nhat_ab_B <- HT (ones_ab_B, pi_B)
			Vhat_Nhat_ab_B <- varest (ones_ab_B, pik = pi_B)
			Vhat_Nhat_ab_A <- varest (ones_ab_A, pik = pi_A)

			for (k in 1:c1) {
	
				data_a_A <- Domains (ysA[,k], domains_A, "a")
				data_ab_A <- Domains (ysA[,k], domains_A, "ab")
				data_b_B <- Domains (ysB[,k], domains_B, "b")
				data_ab_B <- Domains (ysB[,k], domains_B, "ba")

				Yhat_a_A <- HT(data_a_A, pi_A)
				Yhat_ab_A <- HT (data_ab_A, pi_A)
				Yhat_b_B <- HT (data_b_B, pi_B)
				Yhat_ab_B <- HT (data_ab_B, pi_B)

				theta_opt <- (Nhat_a_A * N_B * Vhat_Nhat_ab_B) / (Nhat_a_A * N_B * Vhat_Nhat_ab_B + Nhat_b_B * N_A * Vhat_Nhat_ab_A)
				if (Nhat_a_A == 0 | Nhat_b_B == 0){
		
					SRSPiks_A <- rep(n_A / N_A, n_A)
					SRSPiks_B <- rep(n_B / N_B, n_B)
	
					d_A <- Vhat_Nhat_ab_A / varest (ones_ab_A, pik = SRSPiks_A) 
					d_B <- Vhat_Nhat_ab_B / varest (ones_ab_B, pik = SRSPiks_B)

					eff_n_A <- n_A / d_A
					eff_n_B <- n_B / d_B
		
					theta_opt <- eff_n_A * N_B / (eff_n_A * N_B + eff_n_B * N_A)
				}
				term_a <- theta_opt / N_B + (1 - theta_opt) / N_A
				term_b <- -(1 + theta_opt * Nhat_ab_A / N_B + (1 - theta_opt) * Nhat_ab_B / N_A)
				term_c <- theta_opt * Nhat_ab_A + (1 - theta_opt) * Nhat_ab_B
				Nhat_ab_PML_theta_opt <- (- term_b - sqrt(term_b * term_b - 4 * term_a * term_c)) / (2 * term_a)
				

				muhat_ab <- (n_A / N_A * Yhat_ab_A + n_B / N_B * Yhat_ab_B) / (n_A / N_A * Nhat_ab_A + n_B / N_B * Nhat_ab_B)
				muhat_a <- Yhat_a_A / Nhat_a_A
				muhat_b <- Yhat_b_B / Nhat_b_B

				if (k == 1)
					size_estimation <- sum((N_A - Nhat_ab_PML_theta_opt) * muhat_a, (N_B - Nhat_ab_PML_theta_opt) * muhat_b, Nhat_ab_PML_theta_opt * muhat_ab, rm.na = TRUE)
				else
					total_estimation <- sum((N_A - Nhat_ab_PML_theta_opt) * muhat_a, (N_B - Nhat_ab_PML_theta_opt) * muhat_b, Nhat_ab_PML_theta_opt * muhat_ab, rm.na = TRUE) 
		
				if (k > 1) {
					mean_estimation <- total_estimation / size_estimation
					results[,k-1] <- c(total_estimation, mean_estimation)

					if (!is.null(conf_level)) {
			
						lambdahat <- sum(muhat_ab, - muhat_a, - muhat_b, na.rm = TRUE)
						phihat <- n_A * Nhat_b_B / (n_A * Nhat_b_B + n_B * Nhat_a_A)
						if (is.nan(phihat))
							phihat <- 0
						
						delta_a_A <- Domains (rep (1, n_A), domains_A, "a")
						delta_ab_A <- Domains (rep (1, n_A), domains_A, "ab")
						delta_b_B <- Domains (rep (1, n_B), domains_B, "b")
						delta_ab_B <- Domains (rep (1, n_B), domains_B, "ba")
						
						zA1 <- cbind((ysA[,k] - muhat_a) * delta_a_A, (theta_opt * (ysA[,k] - muhat_ab) + lambdahat * phihat) * delta_ab_A)
						zB1 <- cbind((ysB[,k] - muhat_b) * delta_b_B, ((1 - theta_opt) * (ysB[,k] - muhat_ab) + lambdahat * (1 - phihat)) * delta_ab_B)
						zA <- rowSums(zA1, na.rm = TRUE)
						zB <- rowSums(zB1, na.rm = TRUE)

						Vhat_Yhat_PML <- varest (zA, pik = pi_A) + varest (zB, pik = pi_B)
						total_upper <- total_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_PML)
						total_lower <- total_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(Vhat_Yhat_PML)
						mean_upper <- mean_estimation + qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * Vhat_Yhat_PML)
						mean_lower <- mean_estimation - qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * Vhat_Yhat_PML)
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