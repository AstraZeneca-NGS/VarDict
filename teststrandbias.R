#!/usr/bin/env Rscript

#args <- commandArgs(trailingOnly = TRUE)

#d <- read.table( args[1], sep = "\t", header = F, colClasses=c("character", NA, NA, NA, NA, "character", "character", NA, NA, NA, NA, NA, NA, "character", NA, NA, NA, NA, NA, NA, NA, NA))
d <- tryCatch( {
    read.table( file('stdin'), sep = "\t", header = F, colClasses=c("character", NA, NA, NA, NA, "character", "character", NA, NA, NA, NA, NA, NA, "character", NA, NA, NA, NA, NA, NA, NA, NA))
}, error = function(e) {
    return(NULL)
} )

if (!is.null(d)){
    pvalues <- vector(mode="double", length=dim(d)[1])
    oddratio <- vector(mode="double", length=dim(d)[1])
    
    for( i in 1:dim(d)[1] ) {
        h <- fisher.test(matrix(c(d[i,10], d[i,11], d[i,12], d[i,13]), nrow=2))
        pvalues[i] <- round(h$p.value, 5)
        oddratio[i] <- round(h$estimate, 5)
    }
    write.table(data.frame(d[,1:20], pvalues, oddratio, d[,21:dim(d)[2]]), file = "", quote = F, sep = "\t", eol = "\n", row.names=F, col.names=F)
}
