#!/group/cancer_informatics/tools_resources/R/R-2.15.3/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

d <- read.table( file('stdin'), sep = "\t", header = F, colClasses=c("character", NA, NA, NA, NA, "character", "character", NA, NA, NA, NA, NA, NA, "character", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "character", "character", NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, "character", "character", "character", "character"))

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
    h <- fisher.test(matrix(c(d[i,27], d[i,28], d[i,29], d[i,30]), nrow=2))
    pvalues2[i] <- round(h$p.value, 5)
    oddratio2[i] <- round(h$estimate, 5)
    tref <- if ( d[i,8] - d[i,9] < 0 ) 0 else d[i,8] - d[i,9]
    rref <- if ( d[i,25] - d[i,26] < 0 ) 0 else d[i,25] - d[i,26]
    h <- fisher.test(matrix(c(tref, d[i,9], rref, d[i,26]), nrow=2))
    pvalues[i] <- round(h$p.value, 5)
    oddratio[i] <- round(h$estimate, 5)
}

write.table(data.frame(d[,1:24], pvalues1, oddratio1, d[,25:41], pvalues2, oddratio2, d[, 42:dim(d)[2]], pvalues, oddratio), file = "", quote = F, sep = "\t", eol = "\n", row.names=F, col.names=F)
