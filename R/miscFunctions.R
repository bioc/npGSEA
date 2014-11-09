runNormApprox <- function(xg, y, wg, set){
    betahatG <- .calcBetaHatG (xg, y)
    ThatGw <- .calcThatGw (betahatG, wg)
    varThatGw <- .calcVarThatGw (xg, y, wg)
    stats <- .calcPvaluesNorm (ThatGw, varThatGw)
    return( new("npGSEAResultNorm", 
        geneSetName = setName(set),
        zStat = as.numeric(ThatGw/sqrt(varThatGw)),
        ThatGw = as.numeric(ThatGw),
        varThatGw = as.numeric(varThatGw), 
        pLeft = as.numeric(stats$pLeft), 
        pRight = as.numeric(stats$pRight), 
        pTwoSided = as.numeric(stats$pTwoSided), 
        xSet = xg,
        betaHats = betahatG) )
}

.calcBetaHatG <- function(xg, y){
	nobs <- length(y)
    ###a vector containing each Betahat_g
   	betahat_g <- (xg %*% y)/nobs  
   	return(betahat_g)
}

##test stat for gene set in experiment (norm approx):
.calcThatGw <- function(betahat_g, wg) {
    ##sum of all the Betahat_g's wtd by  their corresponding weights
    ThatGw <- wg %*% betahat_g 
    return(ThatGw)
}

##calculates variance of permuted ThatGw (norm approx)
## based on eq 5, section 3.3 in Larson and Owen 2014
.calcVarThatGw <- function(xg, y, wg) {
    nobs <- length(y)
    XGi <- wg %*% xg 
    XbarGG <- mean(XGi^2)
    mu2 <- mean(y^2)  ##section 3.1
    VarThatGw  <- (mu2/(nobs-1))*XbarGG
    return(VarThatGw)
}


##calculates p-values for our appoximation (norm)
.calcPvaluesNorm <- function(ThatGw, varThatGw){
    Tscaled <- ThatGw/sqrt(varThatGw)
    pL <- pnorm(Tscaled)
    pR <- pnorm(-Tscaled)
    pC <- 2*min(pL, pR)
    return(list(pLeft=pL, pRight=pR, pTwoSided=pC))
}

### calc Beta dist'n based stats
##uses many of the same functions as norm approx
runBetaApprox <- function(xg, y, wg, set, epAdj){    
    #first we get betahats, That and its variance
    betahatG <- .calcBetaHatG (xg, y)
    ThatGw <- .calcThatGw (betahatG, wg)
    varThatGw <- .calcVarThatGw(xg, y, wg)

    ##next we calculate A and B, see section 3.3
    nobs <- length(y)
    XGi <- wg %*% xg 
    sortedXGi <- sort(XGi)
    sortedY <- sort(y)
    sortedYdecreasing <- sort(y, decreasing=TRUE)
    A <- (sortedXGi%*%sortedYdecreasing)/nobs 
    B <- (sortedXGi%*%sortedY)/nobs

    ##next we calculate alpha and beta and our betaStat, see eq. 3
    alpha <- (A/(B-A))*(A*B/varThatGw+1)
    beta <- (-B/(B-A))*(A*B/varThatGw+1)
    betaStat <- (ThatGw-A)/(B-A)
    pvals <- .calcPvaluesBeta (betaStat, alpha, beta, y, epAdj)
    
    return(new("npGSEAResultBeta",
        geneSetName = setName(set),
        betaStat = as.numeric(betaStat),
        ThatGw = as.numeric(ThatGw),
        varThatGw = as.numeric(varThatGw), 
        alpha = as.numeric(alpha), 
        beta = as.numeric(beta), 
        pLeft = as.numeric(pvals$pLeft), 
        pRight = as.numeric(pvals$pRight), 
        pTwoSided = as.numeric(pvals$pTwoSided), 
        xSet = xg,
        betaHats = betahatG)  )
}

##calculates p-values for our appoximation (beta)
##with the epsilon min p-value based on the number of levels of y
.calcPvaluesBeta <- function(betaStat, alpha, beta, y, epAdj){
	epsilon <- 0
	
	if (epAdj==TRUE){
		yFactor <- as.factor(y)
    	yLevels <- levels(yFactor)
   	 	numLevels <- length(yLevels)
    	n <- length(y)
    	epsilon <- 1
    	for (j in 1:numLevels){
    		locLevel <- which(yFactor==yLevels[j])
    		nj <- length(locLevel)
    		epsilon <- epsilon*factorial(nj)
    	}
    	epsilon <- epsilon/factorial(n)
	}
	
	pL <- pbeta(betaStat, alpha, beta)
	pLeft <- epsilon + (1-2*epsilon)*pL
	
    pR <- pbeta(betaStat, alpha, beta, lower.tail=FALSE)
    pRight <- epsilon + (1-2*epsilon)*pR
    
    pC <- 2*min(pLeft, pRight)
    return(list(pLeft=pLeft, pRight=pRight, pTwoSided=pC) )
}




##test stat for gene set in experiment (chisq approx):
runChisqApprox <- function(xg, y, wg, set){
    if( length(y) < 4 ){
        stop("Number of samples is too small 
        for permutation approximation; 
        you need at least 4 samples for the 
        chiSq approximation analysis")
    }
    #Get set statistics
    betahatG <- .calcBetaHatG (xg, y)
    ChatGw <- .calcChatGw(betahatG, wg)
    # Get moment approximation to permutation distribution
    chiSqMoments <- .calcChiSqMoments(xg, y, wg)   
    df <- as.numeric(2*chiSqMoments$mu^2/chiSqMoments$var)
    sigsq <- as.numeric(chiSqMoments$mu/df)
    pQ <- as.numeric(1 - pchisq(ChatGw/sigsq, df) )
   
	##check moments and NaN's
	if (is.nan(pQ)==TRUE){
		print(paste0("The moments for the chisquare approximation (mu= ", chiSqMoments$mu , " and var = ",
		chiSqMoments$var, "), which did not lead to a valid approximation.  This may be a result of a small sample size at one or more levels 
		of the predictor or response variable"))
	}
    
    return(new("npGSEAResultChiSq",
        geneSetName = setName(set),
        chiSqStat = ChatGw/sigsq,
        ChatGw = ChatGw,
        sigmaSq = sigsq, 
        DF= df, 
        pTwoSided = pQ, 
        xSet = xg,
        betaHats = betahatG) )
}

# sum of squared dot products, genes x response y with wts
.calcChatGw <- function(betahat_g, wg) {
    C_Gw <- wg %*%(betahat_g)^2
    return(as.numeric(C_Gw))
}

# Return permutation based moments (chisq approx)
##section 3.1, eq4 of paper
.calcChiSqMoments <- function(xg, y, wg){
    x <- t(xg)
    nobs <- length(y)  ##num of obs
    G = ncol(x)  ##number of genes in G

    ##get A and B matrices
    A <- .getA(nobs)
    B <- .getB()
    AtB = t(A) %*% B  ###2x2 matrix of constants

    ## Precompute X moments
    xgh = xgghh = matrix(NA,G,G)
    for( g in 1:G ){
        for( h in 1:g ){
            xgh[g,h]   = xgh[h,g]   = mean( x[,g] * x[,h] )
            xgghh[g,h] = xgghh[h,g] = mean( x[,g]^2 * x[,h]^2 )
        }
    }

    ## Precompute Y moments
    ymoments <- matrix( 0,2,1 )
    mu2 <- mean(y^2)
    ymoments[1,1] <- mu2^2  
    ymoments[2,1] <- mean(y^4)  #mu4

    ##Calculate mean square (lemma 1)
    meansquare <- rep(0,G)
    for( g in 1:G ){
        meansquare[g] <- mu2 * xgh[g,g]/(nobs-1)
    }

    ##Calculate cov square (lemma 2, collary 2)
    covsquare <- matrix(0,G,G)
    xmoments <- matrix( 0,2,1 )
    for( g in 1:G ){
        for( h in 1:g ){
            xmoments[1,1]  <- (xgh[g,g]*xgh[h,h]+2*xgh[g,h]^2)/nobs^2
            xmoments[2,1]  <- xgghh[g,h]/nobs^3
            covsquare[g,h] = covsquare[h,g] = 
                t(ymoments) %*% AtB %*% xmoments - 
                ymoments[1,1] * xgh[g,g] * xgh[h,h] /(nobs-1)^2
        }
    }
  
    ##Calculate E(C_tildeGw) and var(C_tildeGw) (eq 4, section 3.1)
    mu <- wg%*% meansquare 
    var <- (wg%*% covsquare)%*%wg
    return( list(mu=mu, var = var) )
}

## Calculate matrix A from Lemma 2, section 3.1
.getA <- function(n){
    A = matrix(0, 5, 2)
    A[,1] = c( 0, 0, n/(n-1),-n/((n-1)*(n-2)), 3*n/((n-1)*(n-2)*(n-3)) )
    A[,2] = c(1,-1/(n-1),-1/(n-1),2/((n-1)*(n-2)),-6/((n-1)*(n-2)*(n-3)))
    return(A)
}
## Calculate matrix B from Lemma 2, section 3.1
.getB <- function(){
    B <- matrix(0, 5, 2)
    B[,1] <- c(0, 0, 1, -2, 1)
    B[,2] <- c(1, -4, -3, 12, 6)
    return(B)
}