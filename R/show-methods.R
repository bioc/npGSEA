## show 
setMethod("show", "npGSEAResultNorm", 
    function(object){
        cat(paste0( "Normal Approximation for ", object@geneSetName, "\n"))
        cat(paste0("T_Gw = ", signif(object@ThatGw,3), "\n") )
        cat(paste0("var(T_Gw) = ", signif(object@varThatGw,3), "\n") )
        cat(paste0("pLeft = ", signif(object@pLeft,3), 
                ", pRight = ", signif(object@pRight,3), 
                ", pTwoSided = ",  signif(object@pTwoSided,3), "\n") )
    }
)
setMethod("show", "npGSEAResultBeta", 
    function(object){
        cat(paste0( "Beta Approximation for ", object@geneSetName, "\n"))
        cat(paste0("T_Gw = ", signif(object@ThatGw,3), "\n"))
        cat(paste0("var(T_Gw) = ", signif(object@varThatGw,3), "\n"))
        cat(paste0("pLeft = ", signif(object@pLeft,3), 
                ", pRight = ", signif(object@pRight,3), 
                ", pTwoSided = ",  signif(object@pTwoSided,3), "\n"))
    }
)

setMethod("show", "npGSEAResultChiSq",
    function(object) {
        cat (paste0( "Chi-sq Approximation for ", object@geneSetName), "\n" )
        cat (paste0("C_Gw = ", signif(object@ChatGw,3), "\n"))
        cat (paste0("df = ", signif(object@DF,3), 
            ", sigmaSq = ", signif(object@sigmaSq,3), "\n") )
        cat (paste0("pTwoSided = ",  signif(object@pTwoSided,3), "\n") )
    }
)