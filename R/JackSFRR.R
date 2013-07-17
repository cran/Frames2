#' @name JackSFRR
#' @aliases JackSFRR
#' @title Confidence intervals for raking ratio estimator based on jackknife method
#' 
#' @description Calculates confidence intervals for raking ratio estimator using jackknife procedure
#' 
#' @usage JackSFRR(ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A, 
#' N_B, conf_level, sdA = "srs", sdB = "srs", strA = NULL, strB = NULL, clusA = NULL,   
#' clusB = NULL, NhA = NULL, NhB = NULL, NclushA = NULL, NclushB = NULL, fcpA = FALSE,
#' fcpB = FALSE)
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
#' @param N_A A numeric value indicating the size of frame A
#' @param N_B A numeric value indicating the size of frame B
#' @param conf_level A numeric value indicating the confidence level for the confidence intervals.
#' @param sdA (Optional) A character vector indicating the sampling design considered in frame A. Possible values are "srs" (simple random sampling without replacement), "pps" (probabilities proportional to size sampling), "str" (stratified sampling), "clu" (cluster sampling) and "strclu" (stratified cluster sampling). Default is "srs".
#' @param sdB (Optional) A character vector indicating the sampling design considered in frame B. Possible values are "srs" (simple random sampling without replacement), "pps" (probabilities proportional to size sampling), "str" (stratified sampling), "clu" (cluster sampling) and "strclu" (stratified cluster sampling). Default is "srs".
#' @param strA (Optional) A numeric vector indicating the stratum each unit in frame A belongs to, if a stratified sampling or a stratified cluster sampling has been considered in frame A.
#' @param strB (Optional) A numeric vector indicating the stratum each unit in frame B belongs to, if a stratified sampling or a stratified cluster sampling has been considered in frame B.
#' @param clusA (Optional) A numeric vector indicating the cluster each unit in frame A belongs to, if a cluster sampling or a stratified cluster sampling has been considered in frame A.
#' @param clusB (Optional) A numeric vector indicating the cluster each unit in frame B belongs to, if a cluster sampling or a stratified cluster sampling has been considered in frame B.
#' @param NhA (Optional) A numeric vector indicating the population size in each stratum of frame A. This is only needed when \code{sdA = "str"} and \code{fcpA = TRUE}.
#' @param NhB (Optional) A numeric vector indicating the population size in each stratum of frame B. This is only needed when \code{sdB = "str"} and \code{fcpB = TRUE}.
#' @param NclushA (Optional) A numeric vector indicating the population number of clusters in each stratum of frame A. This is only needed when \code{sdA = "strclu"} and \code{fcpA = TRUE}.
#' @param NclushB (Optional) A numeric vector indicating the population number of clusters in each stratum of frame B. This is only needed when \code{sdB = "strclu"} and \code{fcpB = TRUE}.
#' @param fcpA (Optional) A logic value indicating if a finite population correction factor should be considered in frame A. Default is FALSE.
#' @param fcpB (Optional) A logic value indicating if a finite population correction factor should be considered in frame B. Default is FALSE. 
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
#' @seealso \code{\link{SFRR}}
#' @examples
#' data(DatA)
#' data(DatB) 
#'
#' #Let obtain a 95% jackknife confidence interval for variable Leisure,
#' #supposing a stratified sampling in frame A and a simple random sampling
#' #without replacement in frame B with no finite population correction 
#' #factor in any frame.
#' JackSFRR(DatA$Lei, DatB$Lei, DatA$ProbA, DatB$ProbB, DatA$ProbB, 
#' DatB$ProbA, DatA$Domain, DatB$Domain, 1735, 1191, 0.95, "str", "srs",
#' strA = DatA$Stratum)
#' 
#' #Let check how interval estimation varies when a finite 
#' #population correction factor is considered in both frames.
#' JackSFRR(DatA$Lei, DatB$Lei, DatA$ProbA, DatB$ProbB, DatA$ProbB, 
#' DatB$ProbA, DatA$Domain, DatB$Domain, 1735, 1191, 0.95, "str", "srs", 
#' strA = DatA$Stratum, NhA = c(727, 375, 113, 186, 115, 219), fcpA = TRUE, 
#' fcpB = TRUE)
#' @export
JackSFRR = function (ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A, N_B, conf_level, sdA = "srs", sdB = "srs", strA = NULL, strB = NULL, clusA = NULL, clusB = NULL, NhA = NULL, NhB = NULL, NclushA = NULL, NclushB = NULL, fcpA = FALSE, fcpB = FALSE){

	cnames <- names(ysA)

	ysA <- as.matrix(ysA)
	ysB <- as.matrix(ysB)
	pi_A <- as.matrix(pi_A)
	pi_B <- as.matrix(pi_B)

	c <- ncol(ysA)
	results <- matrix(NA, nrow = 6, ncol = c)
	rownames(results) <- c("Total", "Jack Upper End", "Jack Lower End", "Mean", "Jack Upper End", "Jack Lower End")
	colnames(results) <- cnames

	estimation <- SFRR(ysA, ysB, pi_A, pi_B, pik_ab_B, pik_ba_A, domains_A, domains_B, N_A, N_B)
	size_estimation <- estimation[[2]][1,1] / estimation[[2]][2,1]

	if (sdA == "str"){

		nhA <- table(strA)
		nstrataA <- length(nhA)
		YcstrataA <- matrix(0, nstrataA, c)
		nhA <- c(0, nhA)
		cnhA <- cumsum(nhA)

		for (i in 1:nstrataA){

			k <- 1
			YcA <- matrix(0, nhA[i+1], c)
			for (j in (cnhA[i]+1):cnhA[i+1]){

				if (!is.null(dim(drop(pi_A))))
					YcA[k,] <- SFRR(ysA[-j,], ysB, pi_A[-j,-j], pi_B, pik_ab_B[-j], pik_ba_A, domains_A[-j], domains_B, N_A, N_B)[[2]][1,]
				else
					YcA[k,] <- SFRR(ysA[-j,], ysB, pi_A[-j], pi_B, pik_ab_B[-j], pik_ba_A, domains_A[-j], domains_B, N_A, N_B)[[2]][1,]
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

		if (sdA == "clu"){

			clustersA <- unique(clusA)
			if (!is.null(dim(drop(pi_A))))
				probclustersA <- unique(data.frame(diag(pi_A), clusA))[,1]
			else
				probclustersA <- unique(data.frame(pi_A, clusA))[,1]
			nclustersA <- length(clustersA)
			npsuclustersA <- table(clusA)
			if (any(npsuclustersA < 3))		
				stop("Number of units in any cluster from frame A is less than 3. Variance cannot be computed.")
			
			YcA <- matrix(0, nclustersA, c)

			for (i in 1:nclustersA){
				if (!is.null(dim(drop(pi_A))))
					YcA[i,] <- SFRR(ysA[clusA %in% clustersA[-clustersA[i]],], ysB, pi_A[clusA %in% clustersA[-clustersA[i]],clusA %in% clustersA[-clustersA[i]]], pi_B, pik_ab_B[clusA %in% clustersA[-clustersA[i]]], pik_ba_A, domains_A[clusA %in% clustersA[-clustersA[i]]], domains_B, N_A, N_B)[[2]][1,]
				else
					YcA[i,] <- SFRR(ysA[clusA %in% clustersA[-clustersA[i]],], ysB, pi_A[clusA %in% clustersA[-clustersA[i]]], pi_B, pik_ab_B[clusA %in% clustersA[-clustersA[i]]], pik_ba_A, domains_A[clusA %in% clustersA[-clustersA[i]]], domains_B, N_A, N_B)[[2]][1,]
			}

			YcAMean <- matrix(colMeans(YcA), nclustersA, c, byrow = TRUE)

			if (fcpA == TRUE)
				fA <- 1 - mean(probclustersA)
			else
				fA <- 1

			vjA <- ((nclustersA - 1) / nclustersA) * fA * colSums ((YcA - YcAMean)^2)
		}
		else{
			if (sdA == "strclu"){

				strataA <- unique(strA)
				nstrataA <- length(strataA)
				nhA <- table(strA)
				YcstrataA <- matrix(0, nstrataA, c)
				nhA <- c(0,nhA)
				cnhA <- cumsum(nhA)

				for (i in 1:nstrataA){

					clustersA <- unique(clusA[strA == strataA[i]])
					nclustersA <- length(clustersA)
					clusterpsuA <- table(clusA[strA == strataA[i]])
					if (any(clusterpsuA < 3))		
					stop("Number of units in any cluster from frame A is less than 3. Variance cannot be computed.")
					k <- 1
					YcA <- matrix(0, nclustersA, c)
					for (j in 1:nclustersA){

						if (!is.null(dim(drop(pi_A))))
							YcA[k,] <- SFRR(ysA[clusA %in% clustersA[-clustersA[j]],], ysB, pi_A[clusA %in% clustersA[-clustersA[j]],clusA %in% clustersA[-clustersA[j]]], pi_B, pik_ab_B[clusA %in% clustersA[-clustersA[j]]], pik_ba_A, domains_A[clusA %in% clustersA[-clustersA[j]]], domains_B, N_A, N_B)[[2]][1,]
						else
							YcA[k,] <- SFRR(ysA[clusA %in% clustersA[-clustersA[j]],], ysB, pi_A[clusA %in% clustersA[-clustersA[j]]], pi_B, pik_ab_B[clusA %in% clustersA[-clustersA[j]]], pik_ba_A, domains_A[clusA %in% clustersA[-clustersA[j]]], domains_B, N_A, N_B)[[2]][1,]
						k <- k + 1
					}

					YcAMean <- matrix(colMeans(YcA), nclustersA, c, byrow = TRUE)

					if (fcpA == TRUE){
						fA <- nclustersA/NclushA[i]
					}else {
						fA <- 1
					}
					YcstrataA[i,] <- (nclustersA - 1) / nclustersA * fA * colSums((YcA - YcAMean)^2)
				}
				vjA <- colSums(YcstrataA)
			}else{
				n_A <- nrow(ysA)
				YcA <- matrix(0, n_A, c)
		
				for (i in 1:n_A){
					if (!is.null(dim(drop(pi_A))))
						YcA[i,] <- SFRR(ysA[-i,], ysB, pi_A[-i,-i], pi_B, pik_ab_B[-i], pik_ba_A, domains_A[-i], domains_B, N_A, N_B)[[2]][1,]
					else
						YcA[i,] <- SFRR(ysA[-i,], ysB, pi_A[-i], pi_B, pik_ab_B[-i], pik_ba_A, domains_A[-i], domains_B, N_A, N_B)[[2]][1,]
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
		}
	}
	if (sdB == "str"){

		nhB <- table(strB)
		strataB <- length(nhB)
		YcstrataB <- matrix(0, strataB, c)
		nhB <- c(0,nhB)
		cnhB <- cumsum(nhB)

		for (i in 1:strataB){

			k <- 1
			YcB <- matrix(0, nhB[i+1], c)
			for (j in (cnhB[i]+1):cnhB[i+1]){

				if (!is.null(dim(drop(pi_B))))
					YcB[k,] <- SFRR(ysA, ysB[-j,], pi_A, pi_B[-j,-j], pik_ab_B, pik_ba_A[-j], domains_A, domains_B[-j], N_A, N_B)[[2]][1,]
				else
					YcB[k,] <- SFRR(ysA, ysB[-j,], pi_A, pi_B[-j], pik_ab_B, pik_ba_A[-j], domains_A, domains_B[-j], N_A, N_B)[[2]][1,]
				k <- k + 1
			}

			YcBMean <- matrix(colMeans(YcB), nhB[i+1], c, byrow = TRUE)

			if (fcpB == TRUE)
				fB <- nhB[i+1]/NhB[i]
			else
				fB <- 1
			YcstrataB[i,] <- (nhB[i+1] - 1) / nhB[i+1] * fB * colSums((YcB - YcBMean)^2)
		}
		vjB <- colSums(YcstrataB)	
	}
	else{

		if (sdB == "clu"){

			clustersB <- unique(clusB)
			if (!is.null(dim(drop(pi_B))))
				probclustersB <- unique(data.frame(diag(pi_B), clusB))[,1]
			else
				probclustersB <- unique(data.frame(pi_B, clusB))[,1]
			nclustersB <- length(clustersB)
			npsuclustersB <- table(clusB)
			if (any(npsuclustersB < 3))		
				stop("Number of units in any cluster from frame B is less than 3. Variance cannot be computed.")		

			YcB <- matrix(0, nclustersB, c)

			for (i in 1:nclustersB){
				if (!is.null(dim(drop(pi_B))))
					YcB[i,] <- SFRR(ysA, ysB[clusB %in% clustersB[-clustersB[i]],], pi_A, pi_B[clusB %in% clustersB[-clustersB[i]],clusB %in% clustersB[-clustersB[i]]], pik_ab_B, pik_ba_A[clusB %in% clustersB[-clustersB[i]]], domains_A, domains_B[clusB %in% clustersB[-clustersB[i]]], N_A, N_B)[[2]][1,]
				else
					YcB[i,] <- SFRR(ysA, ysB[clusB %in% clustersB[-clustersB[i]],], pi_A, pi_B[clusB %in% clustersB[-clustersB[i]]], pik_ab_B, pik_ba_A[clusB %in% clustersB[-clustersB[i]]], domains_A, domains_B[clusB %in% clustersB[-clustersB[i]]], N_A, N_B)[[2]][1,]
			}

			YcBMean <- matrix(colMeans(YcB), nclustersB, c, byrow = TRUE)

			if (fcpB == TRUE)
				fB <- 1 - mean(probclustersB)
			else
				fB <- 1

			vjB <- ((nclustersB - 1) / nclustersB) * fB * colSums ((YcB - YcBMean)^2)
		}
		else{

			if(sdB == "strclu"){

				strataB <- unique(strB)
				nstrataB <- length(strataB)
				nhB <- table(strB)
				YcstrataB <- matrix(0, nstrataB, c)
				nhB <- c(0, nhB)
				cnhB <- cumsum(nhB)

				for (i in 1:nstrataB){

					clustersB <- unique(clusB[strB == strataB[i]])
					nclustersB <- length(clustersB)
					clusterpsuB <- table(clusB[strB == strataB[i]])
					if (any(clusterpsuB < 3))		
					stop("Number of units in any cluster from frame B is less than 3. Variance cannot be computed.")
					k <- 1
					YcB <- matrix(0, nclustersB, c)
					for (j in 1:nclustersB){

						if (!is.null(dim(drop(pi_B))))
							YcB[k,] <- SFRR(ysA, ysB[clusB %in% clustersB[-clustersB[j]],], pi_A, pi_B[clusB %in% clustersB[-clustersB[j]],clusB %in% clustersB[-clustersB[j]]], pik_ab_B, pik_ba_A[clusB %in% clustersB[-clustersB[j]]], domains_A, domains_B[clusB %in% clustersB[-clustersB[j]]], N_A, N_B)[[2]][1,]
						else
							YcB[k,] <- SFRR(ysA, ysB[clusB %in% clustersB[-clustersB[j]],], pi_A, pi_B[clusB %in% clustersB[-clustersB[j]]], pik_ab_B, pik_ba_A[clusB %in% clustersB[-clustersB[j]]], domains_A, domains_B[clusB %in% clustersB[-clustersB[j]]], N_A, N_B)[[2]][1,]
						k <- k + 1
					}

					YcBMean <- matrix(colMeans(YcB), nclustersB, c, byrow = TRUE)

					if (fcpB == TRUE){
						fB <- nclustersB/NclushB[i]
					}else {
						fB <- 1
					}
					YcstrataB[i,] <- (nclustersB - 1) / nclustersB * fB * colSums((YcB - YcBMean)^2)
				}
				vjB <- colSums(YcstrataB)
			}else{
				n_B <- nrow(ysB)
				YcB <- matrix(0, n_B, c)

				for (i in 1:n_B){

					if (!is.null(dim(drop(pi_B))))
						YcB[i,] <- SFRR(ysA[-i,], ysB, pi_A[-i,-i], pi_B, pik_ab_B[-i], pik_ba_A, domains_A[-i], domains_B, N_A, N_B)[[2]][1,]
					else
						YcB[i,] <- SFRR(ysA[-i,], ysB, pi_A[-i], pi_B, pik_ab_B[-i], pik_ba_A, domains_A[-i], domains_B, N_A, N_B)[[2]][1,]
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
		}
	}

	VJack_Yhat_SFRR <- vjA + vjB

	results[1,] <- estimation[[2]][1,]
	results[2,] <- estimation[[2]][1,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(VJack_Yhat_SFRR)
	results[3,] <- estimation[[2]][1,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(VJack_Yhat_SFRR)
	results[4,] <- estimation[[2]][2,]
	results[5,] <- estimation[[2]][2,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * VJack_Yhat_SFRR)
	results[6,] <- estimation[[2]][2,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * VJack_Yhat_SFRR)

	return(results)
}