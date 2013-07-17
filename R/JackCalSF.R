#' @name JackCalSF
#' @aliases JackCalSF
#' @title Confidence intervals for single frame calibration estimator based on jackknife method
#' 
#' @description Produces estimates for variance of single frame calibration estimator using Jackknife procedure
#' 
#' @usage JackCalSF(ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, 
#' N_A = NULL, N_B = NULL, N_ab = NULL, XsAFrameA = NULL, XsBFrameA = NULL, 
#' XsAFrameB = NULL, XsBFrameB = NULL, XsT = NULL, XA = NULL, XB = NULL,  
#' X = NULL, met = "linear", conf_level, sdA = "srs", sdB = "srs", nhA = NULL,
#' NhA = NULL, nhB = NULL, NhB = NULL, fcpA = FALSE, fcpB = FALSE)
#' @param ysA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{c} containing information about variable of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{c} containing information about variable of interest from \eqn{s_B}.
#' @param pi_A A numeric vector of length \eqn{n_A} or a square numeric matrix of dimension \eqn{n_A} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param pi_B A numeric vector of length \eqn{n_B} or a square numeric matrix of dimension \eqn{n_B} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param pik_ab_B A numeric vector of size \eqn{n_A} containing first order inclusion probabilities according to sampling desing in frame B for units belonging 
#'  to overlap domain that have been selected in \eqn{s_A}.
#' @param pik_ba_A A numeric vector of size \eqn{n_B} containing first order inclusion probabilities according to sampling desing in frame A for units belonging 
#'  to overlap domain that have been selected in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param N_A (Optional) A numeric value indicating the size of frame A
#' @param N_B (Optional) A numeric value indicating the size of frame B
#' @param N_ab (Optional) A numeric value indicating the size of the overlap domain
#' @param XsAFrameA (Optional) A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_A}.
#' @param XsBFrameA (Optional) A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_B}. For units in domain \eqn{b}, these values are 0.
#' @param XsAFrameB (Optional) A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_A}. For units in domain \eqn{a}, these values are 0.
#' @param XsBFrameB (Optional) A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_B}.
#' @param XsT (Optional) A numeric vector of length \eqn{n} or a numeric matrix or data frame of dimensions \eqn{n} x \eqn{m_T}, with \eqn{m_T} the number of auxiliary variables in both frames, containing auxiliary information for all units in the entire sample \eqn{s = s_A \cup s_B}.
#' @param XA (Optional) A numeric value or vector of length \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, indicating the population totals for the auxiliary variables considered in frame A.
#' @param XB (Optional) A numeric value or vector of length \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, indicating the population totals for the auxiliary variables considered in frame B.
#' @param X (Optional) A numeric value or vector of length \eqn{m_T}, with \eqn{m_T} the number of auxiliary variables in both frames, indicating the population totals for the auxiliary variables considered in both frames.
#' @param met (Optional) A character vector indicating the distance that must be used in calibration process. Possible values are "linear", "raking" and "logit". Default is "linear".
#' @param conf_level A numeric value indicating the confidence level for the confidence intervals.
#' @param sdA (Optional) A character vector indicating the sampling design considered in frame A. Possible values are "srs" (simple random sampling), "pps" (probabilities proportional to size sampling) and "str" (stratified sampling). Default is "srs".
#' @param sdB (Optional) A character vector indicating the sampling design considered in frame B. Possible values are "srs" (simple random sampling), "pps" (probabilities proportional to size sampling) and "str" (stratified sampling). Default is "srs".
#' @param nhA (Optional) A numeric vector indicating the sample size in each stratum, if a stratified sample has been drawn in frame A.
#' @param nhB (Optional) A numeric vector indicating the sample size in each stratum, if a stratified sample has been drawn in frame B.
#' @param NhA (Optional) A numeric vector indicating the population size in each stratum, if a stratified sample has been drawn in frame A.
#' @param NhB (Optional) A numeric vector indicating the population size in each stratum, if a stratified sample has been drawn in frame B.
#' @param fcpA (Optional) A logic value indicating if a finite population correction factor should be considered in frame A. Default is FALSE
#' @param fcpB (Optional) A logic value indicating if a finite population correction factor should be considered in frame B. Default is FALSE 
#' @details Let suppose a non stratified sampling design in frame A and a stratified sampling design in frame B where frame has been divided into L strata and a sample of size \eqn{n_{Bl}} from the \eqn{N_{Bl}} composing the l-th stratum is selected
#' In this context, jackknife variance estimator of a estimator \eqn{\hat{Y}_c} is given by
#'  \deqn{v_J(\hat{Y}_c) = \frac{n_{A}-1}{n_{A}}\sum_{i\in s_A} (\hat{Y}_{c}^{A}(i) -\overline{Y}_{c}^{A})^2 + \sum_{l=1}^{L}\frac{n_{Bl}-1}{n_{Bl}}  \sum_{i\in s_{Bl}} (\hat{Y}_{c}^{B}(lj) -\overline{Y}_{c}^{Bl})^2}
#' with \eqn{\hat{Y}_c^A(i)} the value of estimator \eqn{\hat{Y}_c} after dropping \eqn{i-th} unit from \code{ysA} and \eqn{\overline{Y}_{c}^{A}} the mean of values \eqn{\hat{Y}_c^A(i)}.
#' Similarly, \eqn{\hat{Y}_c^B(lj)} is the value taken by \eqn{\hat{Y}_c} after dropping j-th unit of l-th from sample \code{ysB} and \eqn{\overline{Y}_{c}^{Bl}} is the mean of values \eqn{\hat{Y}_c^B(lj)}.
#' If needed, a finite population correction factor can be included in frames by replacing \eqn{\hat{Y}_{c}^{A}(i)} or \eqn{\hat{Y}_{c}^{B}(lj)} with \eqn{\hat{Y}_{c}^{A*}(i)= \hat{Y}_{c}+\sqrt{1-\overline{\pi}_A} (\hat{Y}_{c}^{A}(i) -\hat{Y}_{c})} or
#' \eqn{\hat{Y}_{c}^{B*}(lj)= \hat{Y}_{c}+\sqrt{1-\overline{\pi}_B} (\hat{Y}_{c}^{B}(lj) -\hat{Y}_{c})}, where \eqn{\overline{\pi}_A = \sum_{i \in s_A}\pi_{iA}/n_A} and \eqn{\overline{\pi}_B = \sum_{j \in s_B}\pi_{jB}/n_B}
#' A confidence interval for any parameter of interest, \eqn{Y} can be calculated, then, using the pivotal method.
#' @return A numeric matrix containing estimations of population total and population mean and their corresponding confidence intervals obtained through jackknife method.
#' @references Wolter, K. M. (2007)
#'  \emph{Introduction to Variance Estimation.}
#'  2nd Edition. Springer, Inc., New York.
#' @seealso \code{\link{CalSF}}
#' @examples
#' data(HouseholdsA)
#' dataA <- attach(HouseholdsA)
#' detach(HouseholdsA)
#' data(HouseholdsB)
#' dataB <- attach(HouseholdsB)
#' detach(HouseholdsB)
#' 
#' #Let obtain a 95% jackknife confidence interval for variable Clothing,
#' #with frame sizes and overlap domain size known, supposing a stratified
#' #sampling in frame A and a simple random sampling in frame B with no
#' #finite population correction factor in any frame
#' JackCalSF(dataA$Clothing, dataB$Clothing, dataA$ProbA, dataB$ProbB, 
#' dataA$ProbB, dataB$ProbA, dataA$Domain, dataB$Domain, N_A = 1735, 
#' N_B = 1191, N_ab = 601, conf_level = 0.95, sdA = "str", sdB = "srs", 
#' nhA = c(15, 20, 15, 20, 15, 20))
#' @export
JackCalSF = function (ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A = NULL, N_B = NULL, N_ab = NULL, XsAFrameA = NULL, XsBFrameA = NULL, XsAFrameB = NULL, XsBFrameB = NULL, XsT = NULL, XA = NULL, XB = NULL, X = NULL, met = "linear", conf_level, sdA = "srs", sdB = "srs", nhA = NULL, NhA = NULL, nhB = NULL, NhB = NULL, fcpA = FALSE, fcpB = FALSE){

	cnames <- names(ysA)

	ysA <- as.matrix(ysA)
	ysB <- as.matrix(ysB)
	pi_A <- as.matrix(pi_A)
	pi_B <- as.matrix(pi_B)
	if (!is.null(XsAFrameA)){
		XsAFrameA <- as.matrix(XsAFrameA)
		XsAFrameB <- as.matrix(XsAFrameB)
	}
	if (!is.null(XsBFrameA)){
		XsBFrameA <- as.matrix(XsBFrameA)
		XsBFrameB <- as.matrix(XsBFrameB)
	}
	if (!is.null(XsT))
		XsT <- as.matrix(XsT)
	
	n_A <- nrow(ysA)
	n_B <- nrow(ysB)
	c <- ncol(ysA)
	results <- matrix(NA, nrow = 6, ncol = c)
	rownames(results) <- c("Total", "Jack Upper End", "Jack Lower End", "Mean", "Jack Upper End", "Jack Lower End")
	colnames(results) <- cnames

	estimation <- CalSF(ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A, N_B, N_ab, XsAFrameA, XsBFrameA, XsAFrameB, XsBFrameB, XsT, XA, XB, X, met)
	size_estimation <- estimation[1,1] / estimation[2,1]

	if (sdA == "str"){

		strataA <- length(nhA)
		YcstrataA <- matrix(0, strataA, c)
		nhA <- c(0,nhA)
		cnhA <- cumsum(nhA)

		for (i in 1:strataA){

			k <- 1
			YcA <- matrix(0, nhA[i+1], c)
			for (j in (cnhA[i]+1):cnhA[i+1]){

				if (!is.null(dim(drop(pi_A))))
					YcA[k,] <- CalSF(ysA[-j,], ysB, pi_A[-j,-j], pi_B, pik_ab_B[-j], pik_ba_A, domains_A[-j], domains_B, N_A, N_B, N_ab, XsAFrameA[-j,], XsBFrameA, XsAFrameB[-j,], XsBFrameB, XsT[-j,], XA, XB, X, met)[1,]
				else
					YcA[k,] <- CalSF(ysA[-j,], ysB, pi_A[-j], pi_B, pik_ab_B[-j], pik_ba_A, domains_A[-j], domains_B, N_A, N_B, N_ab, XsAFrameA[-j,], XsBFrameA, XsAFrameB[-j,], XsBFrameB, XsT[-j,], XA, XB, X, met)[1,]
				k <- k + 1
			}

			YcAMean <- matrix(colMeans(YcA), nhA[i+1], c, byrow = TRUE)

			if (fcpA == TRUE){
				fA <- nhA[i+1]/NhA[i]
			}else {
				fA <- 1
			}
			YcstrataA[i,] <- (nhA[i+1] - 1) / nhA[i+1] * fA * colSums((YcA - YcAMean)^2)
		}
		vjA <- colSums(YcstrataA)	
	}
	else {

		YcA <- matrix(0, n_A, c)
		
		for (i in 1:n_A){

			if (!is.null(dim(drop(pi_A))))
				YcA[i,] <- CalSF(ysA[-i,], ysB, pi_A[-i,-i], pi_B, pik_ab_B[-i], pik_ba_A, domains_A[-i], domains_B, N_A, N_B, N_ab, XsAFrameA[-i,], XsBFrameA, XsAFrameB[-i,], XsBFrameB, XsT[-i,], XA, XB, X, met)[1,]
			else
				YcA[i,] <- CalSF(ysA[-i,], ysB, pi_A[-i], pi_B, pik_ab_B[-i], pik_ba_A, domains_A[-i], domains_B, N_A, N_B, N_ab, XsAFrameA[-i,], XsBFrameA, XsAFrameB[-i,], XsBFrameB, XsT[-i,], XA, XB, X, met)[1,]		
		}

		YcAMean <- matrix(colMeans(YcA), n_A, c, byrow = TRUE)
	
		if (fcpA == TRUE)
			if (!is.null(dim(drop(pi_A))))
				fA <- 1 - mean(diag(pi_A))
			else
				fA <- 1 - mean(pi_A)
		else
			fA <- 1

		vjA <- ((n_A - 1) / n_A) * fA * colSums ((YcA - YcAMean)^2)	
	}

	if (sdB == "str"){

		strataB <- length(nhB)
		YcstrataB <- matrix(0, strataB, c)
		nhB <- c(0,nhB)
		cnhB <- cumsum(nhB)

		for (i in 1:strataB){

			k <- 1
			YcB <- matrix(0, nhB[i+1], c, byrow = TRUE)
			for (j in (cnhB[i]+1):cnhB[i+1]){

				if (!is.null(dim(drop(pi_B))))
					YcB[k,] <- CalSF(ysA, ysB[-j,], pi_A, pi_B[-j,-j], pik_ab_B, pik_ba_A[-j], domains_A, domains_B[-j], N_A, N_B, N_ab, XsAFrameA, XsBFrameA[-j,], XsAFrameB, XsBFrameB[-j,], XsT[-j,], XA, XB, X, met)[1,]		
				else 
					YcB[k,] <- CalSF(ysA, ysB[-j,], pi_A, pi_B[-j], pik_ab_B, pik_ba_A[-j], domains_A, domains_B[-j], N_A, N_B, N_ab, XsAFrameA, XsBFrameA[-j,], XsAFrameB, XsBFrameB[-j,], XsT[-j,], XA, XB, X, met)[1,]
				k <- k + 1
			}

			YcBMean <- matrix(colMeans(YcB), nhB[i+1], c)

			if (fcpB == TRUE)
				fB <- nhB[i+1]/NhB[i]
			else
				fB <- 1
			YcstrataB[i,] <- (nhB[i+1] - 1) / nhB[i+1] * fB * colSums((YcB - YcBMean)^2)
		}
		vjB <- colSums(YcstrataB)	
	}
	else{

		YcB <- matrix(0, n_B, c)

		for (i in 1:n_B){

	        	if (!is.null(dim(drop(pi_B))))
				YcB[i,] <- CalSF(ysA, ysB[-i,], pi_A, pi_B[-i,-i], pik_ab_B, pik_ba_A[-i], domains_A, domains_B[-i], N_A, N_B, N_ab, XsAFrameA, XsBFrameA[-i,], XsAFrameB, XsBFrameB[-i,], XsT[-i,], XA, XB, X, met)[1,]
			else
				YcB[i,] <- CalSF(ysA, ysB[-i,], pi_A, pi_B[-i], pik_ab_B, pik_ba_A[-i], domains_A, domains_B[-i], N_A, N_B, N_ab, XsAFrameA, XsBFrameA[-i,], XsAFrameB, XsBFrameB[-i,], XsT[-i,], XA, XB, X, met)[1,]
		}

		YcBMean <- matrix(colMeans(YcB), n_B, c, byrow = TRUE)

		if (fcpB == TRUE)
			if (!is.null(dim(drop(pi_B))))
				fB <- 1-mean(diag(pi_B))
			else
				fB <- 1-mean(pi_B)
		else
			fB <- 1

		vjB <- ((n_B - 1) / n_B) * fB * colSums((YcB - YcBMean)^2)
	}

	VJack_Yhat_CalSF <- vjA + vjB

	results[1,] <- estimation[1,]
	results[2,] <- estimation[1,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(VJack_Yhat_CalSF)
	results[3,] <- estimation[1,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(VJack_Yhat_CalSF)
	results[4,] <- estimation[2,]
	results[5,] <- estimation[2,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * VJack_Yhat_CalSF)
	results[6,] <- estimation[2,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * VJack_Yhat_CalSF)
	
	return(results)
}