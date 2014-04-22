setMethod("geneSetName", "npGSEAResultBeta", 
    function(object){object@geneSetName}
)
setMethod("geneSetName", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@geneSetName)}
)

setMethod("stat", "npGSEAResultBeta", 
    function(object){object@ThatGw}
)
setMethod("stat", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@ThatGw)}
)

setMethod("sigmaSq", "npGSEAResultBeta", 
    function(object){object@varThatGw}
)
setMethod("sigmaSq", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@sigmaSq)}
)

setMethod("pTwoSided", "npGSEAResultBeta", 
    function(object){object@pTwoSided}
)
setMethod("pTwoSided", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@pTwoSided)}
)

setMethod("pLeft", "npGSEAResultBeta", 
    function(object){object@pLeft}
)
setMethod("pLeft", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@pLeft)}
)

setMethod("pRight", "npGSEAResultBeta", 
    function(object){object@pRight}
)
setMethod("pRight", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@pRight)}
)

setMethod("xSet", "npGSEAResultBeta", 
    function(object){object@xSet}
)
setMethod("xSet", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@xSet)}
)

setMethod("betaStat", "npGSEAResultBeta", 
    function(object){object@betaStat}
)
setMethod("betaStat", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@betaStat)}
)

setMethod("alphaValue", "npGSEAResultBeta", 
    function(object){object@alpha}
)
setMethod("alphaValue", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@alpha)}
)

setMethod("betaValue", "npGSEAResultBeta", 
    function(object){object@beta}
)
setMethod("betaValue", "npGSEAResultBetaCollection", 
    function(object){lapply(object, function(y) y@beta)}
)