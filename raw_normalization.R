setwd("D:\\STUDY\\szbl\\metastatic cell genomic\\code\\Primary_CNS")


# Raw data extraction
library(GEOquery)
gse14682 <- getGEO("GSE14682", getGPL = TRUE, GSEMatrix = FALSE)
eset14682 <- getGEO("GSE14682")
raw14682 = exprs(eset14682[[1]])
head(raw14682)

gse14683 <- getGEO("GSE14683", getGPL = TRUE, GSEMatrix = FALSE)
eset14683 <- getGEO("GSE14683")
raw14683 = exprs(eset14683[[1]])
head(raw14683)


for (i in 1:60) {
  gsm <- Table(GSMList(gse14682)[[i]])
  gsm[,"raw"] = as.numeric(gsub(",","",gsm$raw))
  gsm = gsm[match(rownames(raw14682),gsm$ID_REF),]
  raw14682[,names(GSMList(gse14682))[i]] = as.numeric(gsm[,"raw"])
}
head(raw14682)


for (i in 1:length(names(GSMList(gse14683)))) {
  gsm <- Table(GSMList(gse14683)[[i]])
  gsm[,"raw"] = as.numeric(gsub(",","",gsm$raw))
  gsm = gsm[match(rownames(raw14683),gsm$ID_REF),]
  raw14683[,names(GSMList(gse14683))[i]] = as.numeric(gsm[,"raw"])
}
head(raw14683)

# Normalization
library(edgeR)
dge14682 <- DGEList(counts=raw14682)
dge14682 <- calcNormFactors(dge14682)
exp14682 <- cpm(dge14682, log = T, prior.count = 1)
exp14682[1:4, 1:4]
boxplot(exp14682)

dge14683 <- DGEList(counts=raw14683)
dge14683 <- calcNormFactors(dge14683)
exp14683 <- cpm(dge14683, log = T, prior.count = 1)
exp14683[1:4, 1:4]
boxplot(exp14683)

save(exp14682,exp14683,file = "normalized_exp.Rdata")









