setMethod("geneSetName", "npGSEAResultChiSq", 
    function(object){object@geneSetName}
)
setMethod("geneSetName", "npGSEAResultChiSqCollection", 
    function(object){lapply(object, function(y) y@geneSetName)}
)

setMethod("stat", "npGSEAResultChiSq", 
    function(object){object@ChatGw}
)
setMethod("stat", "npGSEAResultChiSqCollection", 
    function(object){lapply(object, function(y) y@ChatGw)}
)

setMethod("sigmaSq", "npGSEAResultChiSq", 
    function(object){object@sigmaSq}
)
setMethod("sigmaSq", "npGSEAResultChiSqCollection", 
    function(object){lapply(object, function(y) y@sigmaSq)}
)

setMethod("pTwoSided", "npGSEAResultChiSq", 
    function(object){object@pTwoSided}
)
setMethod("pTwoSided", "npGSEAResultChiSqCollection", 
    function(object){lapply(object, function(y) y@pTwoSided)}
)

setMethod("xSet", "npGSEAResultChiSq", 
    function(object){object@xSet}
)
setMethod("xSet", "npGSEAResultChiSqCollection", 
    function(object){lapply(object, function(y) y@xSet)}
)

setMethod("chiSqStat", "npGSEAResultChiSq", 
    function(object){object@chiSqStat}
)
setMethod("chiSqStat", "npGSEAResultChiSqCollection", 
    function(object){lapply(object, function(y) y@chiSqStat)}
)

setMethod("DF", "npGSEAResultChiSq", 
    function(object){object@DF}
)
setMethod("DF", "npGSEAResultChiSqCollection", 
    function(object){lapply(object, function(y) y@DF)}
)
