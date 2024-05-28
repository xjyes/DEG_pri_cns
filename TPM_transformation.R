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


# Select useful datasets
load("compare.Rdata")
validate = ifelse(compare$Dataset_id == "GSE14682",
                  as.character(compare$Sample_id),"")
validate = subset(validate,validate != "")
validate
raw14682 = raw14682[,c(validate)]
raw14682 <- as.data.frame(raw14682)

validate = ifelse(compare$Dataset_id == "GSE14683",
                  as.character(compare$Sample_id),"")
validate = subset(validate,validate != "")
validate
raw14683 = raw14683[,c(validate)]
raw14683 = as.data.frame(raw14683)

# Convert to gene name
gpl8128 <- getGEO('GPL8128', destdir = ".")

id8128 <- gpl8128@dataTable@table[,c("ID", "SYMBOL")]


raw14682 = raw14682[rownames(raw14682) %in% id8128$ID,]
id8128 = id8128[match(rownames(raw14682), id8128$ID),]
tmp = by(raw14682,
         id8128$SYMBOL,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
raw14682 = raw14682[rownames(raw14682) %in% probes, ]
rownames(raw14682) = id8128[match(rownames(raw14682),id8128$ID),2]

raw14683 = raw14683[rownames(raw14683) %in% id8128$ID,]
id8128 = id8128[match(rownames(raw14683), id8128$ID),]
tmp = by(raw14683,
         id8128$SYMBOL,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
raw14683 = raw14683[rownames(raw14683) %in% probes, ]
rownames(raw14683) = id8128[match(rownames(raw14683),id8128$ID),2]

raw63668 = raw63668[rownames(raw63668) %in% id16686$ID,]
head(raw63668)
id16686 = id16686[match(rownames(raw63668), id16686$ID),]
tmp = by(raw63668,
         id16686$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
raw63668 = raw63668[rownames(raw63668) %in% probes, ]
rownames(raw63668) = id16686[match(rownames(raw63668),id16686$ID),2]

# Merge gene sets
library("dplyr")
raw14682$symbol = rownames(raw14682)
raw14683$symbol = rownames(raw14683)
merge_eset = inner_join(raw14682,raw14683, by = "symbol")
dim(merge_eset)
rownames(merge_eset) = merge_eset$symbol
merge_eset = merge_eset[, -grep("symbol",colnames(merge_eset))]

save(merge_eset, file = "merge_eset.Rdata")

# TPM transformation
rm(list = ls())
load("merge_eset.Rdata")
library(GenomicFeatures)
txdb <- makeTxDbFromGFF("Homo_sapiens.GRCh38.107.gtf", format = "gtf")
exons_gene <- exonsBy(txdb, by = "gene")
exons_gene_lens <- sum(width(reduce(exons_gene)))
exons_gene_lens[1:10]

gene_length <- sapply(exons_gene_lens, function(x){x})
id_length <- as.data.frame(gene_length)
id_length$gene_ID <- rownames(id_length)
id_length <- as.data.frame(id_length[sort(rownames(id_length)),])
head(id_length)

GN <- read.table("mart_export.txt", head = T, fill = T)
head(GN)
rownames(GN) <- GN$Gene_stable_ID
GN <- GN[sort(rownames(GN)),]

GN = GN[rownames(GN) %in% rownames(id_length),]
id_GN_length <- cbind(GN, id_length)
id_GN_length = id_GN_length[, -grep("Gene_stable_ID",colnames(id_GN_length))]
id_GN_length = id_GN_length[, -grep("gene_id", colnames(id_GN_length))]
colnames(id_GN_length) <- c( "gene_name", "length")
head(id_GN_length)

length <- id_GN_length
rawcount <- merge_eset
rawcount$gene_name <- rownames(rawcount)

library(dplyr)
# length <- distinct(rawcount, gene_name, .keep_all = T)
mergecount <- merge(rawcount, length, by = "gene_name")
FPKMlength <- mergecount[,c(1, ncol(mergecount))]
FPKMcount <- mergecount[, -c(ncol(mergecount))]
rownames(FPKMcount) <- FPKMcount$gene_name
FPKMcount <- FPKMcount[,-1]

# Count TPM
kb <- FPKMlength$length /1000
head(kb)
rpk <- FPKMcount / kb
head(rpk)
tpm <- t(t(rpk)/colSums(rpk) * 1000000)
head(tpm)
dftpm <- as.data.frame(tpm)
write.csv(dftpm, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\immune filtration\\TPM.csv")

