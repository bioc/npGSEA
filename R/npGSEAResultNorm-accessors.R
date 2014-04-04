setMethod("geneSetName", "npGSEAResultNorm", 
    function(object){object@geneSetName}
)
setMethod("geneSetName", "npGSEAResultNormCollection", 
    function(object){lapply(object, function(y) y@geneSetName)}
)

setMethod("stat", "npGSEAResultNorm", 
    function(object){object@ThatGw}
)
setMethod("stat", "npGSEAResultNormCollection", 
    function(object){lapply(object, function(y) y@ThatGw)}
)

setMethod("sigmaSq", "npGSEAResultNorm", 
    function(object){object@varThatGw}
)
setMethod("sigmaSq", "npGSEAResultNormCollection", 
    function(object){lapply(object, function(y) y@varThatGw)}
)

setMethod("pTwoSided", "npGSEAResultNorm", 
    function(object){object@pTwoSided}
)
setMethod("pTwoSided", "npGSEAResultNormCollection", 
    function(object){lapply(object, function(y) y@pTwoSided)}
)

setMethod("pLeft", "npGSEAResultNorm", 
    function(object){object@pLeft}
)
setMethod("pLeft", "npGSEAResultNormCollection", 
    function(object){lapply(object, function(y) y@pLeft)}
)

setMethod("pRight", "npGSEAResultNorm", 
    function(object){object@pRight}
)
setMethod("pRight", "npGSEAResultNormCollection", 
    function(object){lapply(object, function(y) y@pRight)}
)

setMethod("xSet", "npGSEAResultNorm", 
    function(object){object@xSet}
)
setMethod("xSet", "npGSEAResultNormCollection", 
    function(object){lapply(object, function(y) y@xSet)}
)

setMethod("zStat", "npGSEAResultNorm", 
    function(object){object@zStat}
)
setMethod("zStat", "npGSEAResultNormCollection", 
    function(object){lapply(object, function(y) y@zStat)}
)







