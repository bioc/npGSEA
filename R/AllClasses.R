setClass("npGSEAResultNorm",
    representation(
    geneSetName = "character",
    zStat = "numeric",
    ThatGw = "numeric",
    varThatGw = "numeric",
    pLeft="numeric",
    pRight="numeric",
    pTwoSided="numeric",
    xSet = "matrix",
    betaHats = "matrix")
)

setClass("npGSEAResultNormCollection",
    contains="list",
    validity = function(object) {
        msg <- NULL
        if (!all(sapply(object, is, "npGSEAResultNorm")))
            msg <- c(msg, 
                "members must all be 'npGSEAResultNorm' classes")
        if (!is.null(msg))
            msg
        else
            TRUE
    }
)

setClass("npGSEAResultBeta",
    representation(
    geneSetName = "character",
    betaStat="numeric",
    ThatGw = "numeric",
    varThatGw = "numeric",
    alpha="numeric",
    beta="numeric",
    pLeft="numeric",
    pRight="numeric",
    pTwoSided="numeric",
    xSet="matrix",
    betaHats = "matrix")
)

setClass("npGSEAResultBetaCollection",
    contains="list",
    validity = function(object) {
        msg <- NULL
        if (!all(sapply(object, is, "npGSEAResultBeta")))
            msg <- c(msg, 
                "members must all be 'npGSEAResultBeta' classes")
        if (!is.null(msg))
            msg
        else
            TRUE
         }
)

setClass("npGSEAResultChiSq",
    representation(
    geneSetName = "character",
    chiSqStat="numeric",
    ChatGw="numeric",
    sigmaSq="numeric",
    DF="numeric",
    pTwoSided="numeric",
    xSet="matrix",
    betaHats = "matrix")
)

setClass("npGSEAResultChiSqCollection",
    contains="list",
    validity = function(object) {
        msg <- NULL
        if (!all(sapply(object, is, "npGSEAResultChiSq")))
            msg <- c(msg, 
                "members must all be 'npGSEAResultChiSq' classes")
        if (!is.null(msg))
            msg
        else
            TRUE
    }
)

