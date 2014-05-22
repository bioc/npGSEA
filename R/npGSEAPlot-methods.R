setMethod("npGSEAPlot", "npGSEAResultNorm", 
    function(object){
        xrange <- c(-4,4)
        obsZStat <- object@zStat
        if(obsZStat >= 4) {xrange <- c( -(obsZStat+1), (obsZStat+1) ) } 
        if(obsZStat <= -4) {xrange <- c( (obsZStat-1), -(obsZStat-1) ) }
        plot(function(x) dnorm(x),  
            main = paste0("Standard Normal Distribution"), 
            xlim=xrange, xlab="Z values", ylab="Density")
        abline(v = obsZStat, col="red")
    }
)

setMethod("npGSEAPlot", "npGSEAResultBeta", 
    function(object){
        xrange <- c(0,1)
        obsBStat <- object@betaStat
        alph <- object@alpha
        bet <- object@beta
        if(obsBStat >= 1) { xrange<- c( -(obsBStat+1),(obsBStat+1) ) }
        if(obsBStat <= 0) { xrange<- c( (obsBStat-1),-(obsBStat-1) ) }
        plot(function(x) dbeta(x, alph, bet), 
                main = paste0("Beta distribution with \n alpha = ", 
                    round(alph, 2), " and beta = ", round(bet,2)), 
                xlim = xrange, xlab = "Beta values", ylab = "Density")
        abline(v = obsBStat, col="red")
    }
)
setMethod("npGSEAPlot", "npGSEAResultChiSq", 
    function(object){
        xrange <- c(0,4)
        obsStat <- object@chiSqStat
        df <- object@DF
        if(obsStat >= 4) {xrange <- c(0, (obsStat+1) ) }
        plot(function(x) dchisq(x, df),  
            main = paste0("Chi-sq distribution with \n df=", 
                round(df, 2)), xlim=xrange, xlab="Chi-sq values", ylab="Density")
        abline(v = obsStat, col="red")
    }
)



