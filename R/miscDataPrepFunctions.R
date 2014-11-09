.prepXYZ <- function(x, y,  z=NULL, set) {
    ##get locale of set genes in experiment
    inset <- .locGenes(x, set)
    ##select out genes in set
    xset <- x[which(inset==1),]  ##rows are genes
    ##scale x to sum to zero
    xg <- .center(xset)
    
    .checkxg(xg)

    ##adjust for covariates for each gene and sample, if necessary
    ##z is our covars
   if(is.null(z) == FALSE ){
        gnames <- rownames(xg)
        xg <- t(apply(xg, 1, .getResids, z))
        rownames(xg) <- gnames
        y <- .getResids(y, z)
    }   
    return(list(xg=xg, y=y, inset=inset))
}

##check xg so that there is an error if 
##X_G takes only one value, 
##or if X_G takes only two values but one or both of those values appears only 1 time
.checkxg <- function(x){
	XGi <- apply(x, 2, sum)  ##sum the col'ns
	uniqueXG <- unique(XGi)	
    numVals = length(uniqueXG)
    if (numVals == 1){
    	stop( "You need at least two different values (across all samples) for the sum of the Xg's (for all g in G).") 
    }
    if (numVals == 2){
    	locVal1 = which(XGi==uniqueXG[1])
    	if (length(locVal1)==1){
    		stop( "If the sum of the Xg's (for all g in G) takes only two distinct values, then one of these values needs to appear more than once.")
    	}
    	locVal2 = which(XGi==uniqueXG[2])
    	if (length(locVal2)==1){
    		stop( "If the sum of the Xg's (for all g in G) takes only two distinct values, then one of these values needs to appear more than once.")
    	}
    }
}


##take into account any weights, or assign them all to 1, if none given
.prepW <- function(w, set, inset){
    if( (is.null(w) == FALSE) && ( length(w)!=sum(inset) ) ){
        stop( paste0 ("You must have the same number of weights as genes 
            in your set in this experiment.  ", 
            setName(set), " does not have the same number of 
            genes as there are 
            weights.\nRun getIncidence(rownames(x), set) 
            to determine which 
            genes are in this set and experiment." ) ) 
        }
    if (is.null(w) )  { w <- rep(1, sum(inset)) }
    return(w)
}

##get locale of set genes in experiment
##x is our x matrix (rows are genes, columns are samples)
##set is a GeneSet
.locGenes <- function(x, set) {
    incidence <- getIncidence(rownames(x), set)
    inset <- incidence$inSet
    if(sum(inset) == 0)
        {stop( paste0("No genes in ", setName(set), 
            " are in this experiment.  Make sure that the row names of x 
            are the same type as geneIds(set)") )
        }
    if(sum(inset) == 1)
        {stop( paste0("Only one gene in your set ", setName(set), 
            " is in this experiment.  Please choose a larger 
            set for this analysis"))
       }    
    return(inset)
}
        


##want sum(y)=0
##and want to make sure that there are at least two obs in each level of y
.adjustY <- function(y) {
	yFactor = as.factor(y)
    yLevels = levels(yFactor)
    numLevels=length(yLevels)
    if (numLevels < 3){
   	 	for (l in 1:numLevels){
    		locLevel = which(yFactor==yLevels[l])
    		if (length(locLevel)<2){
    			stop( paste0("Fewer than two observations in group ", yLevels[l], 
            " are in this experiment.  Make sure that there are at least two 
            observations in each treatment group if you only have two treatment groups.") ) 
    		}
    	}
	}
          
    y <- as.numeric(y)
    ynew <- y-mean(y)
    return(as.numeric(ynew))
}

# center rows of x around their mean
.center <- function(x) {
  xscaled <- scale(t(x), center=TRUE, scale=FALSE)
  return(t(xscaled))
}

##adjusts for covariates by taking residuals of fit
.getResids<-function(y1, x1) {
    return(resid(lm(y1~x1)))
}

