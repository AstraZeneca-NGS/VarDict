#!/usr/bin/env Rscript

#args <- commandArgs(trailingOnly = TRUE)

d <- tryCatch( {
    d <- read.table( file('stdin'), sep = "\t", header = F, colClasses=c("character", NA, NA, NA, NA, "character", "character", NA, NA, NA, NA, NA, NA, "character", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "character", NA, "character",  NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "character", "character", "character", "character"))
}, error = function(e) {
    return(NULL)
} )


if (!is.null(d)){
    pvalues1 <- vector(mode="double", length=dim(d)[1])
    oddratio1 <- vector(mode="double", length=dim(d)[1])
    pvalues2 <- vector(mode="double", length=dim(d)[1])
    oddratio2 <- vector(mode="double", length=dim(d)[1])
    pvalues <- vector(mode="double", length=dim(d)[1])
    oddratio <- vector(mode="double", length=dim(d)[1])

    for( i in 1:dim(d)[1] ) {
	h <- fisher.test(matrix(c(d[i,10], d[i,11], d[i,12], d[i,13]), nrow=2))
	pvalues1[i] <- round(h$p.value, 5)
	oddratio1[i] <- round(h$estimate, 5)
	h <- fisher.test(matrix(c(d[i,28], d[i,29], d[i,30], d[i,31]), nrow=2))
	pvalues2[i] <- round(h$p.value, 5)
	oddratio2[i] <- round(h$estimate, 5)
	tref <- if ( d[i,8] - d[i,9] < 0 ) 0 else d[i,8] - d[i,9]
	rref <- if ( d[i,26] - d[i,27] < 0 ) 0 else d[i,26] - d[i,27]
	h <- fisher.test(matrix(c(tref, d[i,9], rref, d[i,27]), nrow=2))
	pvalues[i] <- round(h$p.value, 5)
	oddratio[i] <- round(h$estimate, 5)
    }

    write.table(data.frame(d[,1:25], pvalues1, oddratio1, d[,26:43], pvalues2, oddratio2, d[, 44:dim(d)[2]], pvalues, oddratio), file = "", quote = F, sep = "\t", eol = "\n", row.names=F, col.names=F)
}
