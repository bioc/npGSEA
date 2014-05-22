setMethod("pValues", "npGSEAResultNorm", 
    function(object){
        cat(paste0("pLeft = ", signif(object@pLeft,3), 
                ", pRight = ", signif(object@pRight,3), 
                ", pTwoSided = ",  signif(object@pTwoSided,3), "\n") )
    }
)

setMethod("pValues", "npGSEAResultBeta", 
    function(object){
        cat(paste0("pLeft = ", signif(object@pLeft,3), 
                ", pRight = ", signif(object@pRight,3), 
                ", pTwoSided = ",  signif(object@pTwoSided,3), "\n"))
    }
)

setMethod("pValues", "npGSEAResultChiSq",
    function(object) {
        cat (paste0("pTwoSided = ",  signif(object@pTwoSided,3), "\n") )
    }
)