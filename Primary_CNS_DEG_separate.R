setwd("/Users/xujingyi/Documents/szbl/metastatic cell genomic/code/Primary_CNS")

load("normalized_exp.Rdata")
data3521 <- read.table(file = "/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/Primary_CNS/RNA-seq/expression matrix file/GSE3521-GPL1390_series_matrix.txt",
                       header = T, sep = "\t", quote = "", fill = T, comment.char = "!")
data10893 <- read.table(file = "/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/Primary_CNS/RNA-seq/expression matrix file/GSE10893-GPL1390_series_matrix.txt",
                        header = T, sep = "\t", quote = "", fill = T, comment.char = "!")
data14682 <- exp14682
data14683 <- exp14683
data37407 <- read.table(file = "/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/Primary_CNS/RNA-seq/expression matrix file/GSE37407_series_matrix.txt",
                        header = T, sep = "\t", quote = "", fill = T, comment.char = "!")
data63668 <- read.table(file = "/Users/xujingyi/Documents/szbl/metastatic cell genomic/dataset/Primary_CNS/RNA-seq/expression matrix file/GSE63668_series_matrix.txt",
                        header = T, sep = "\t", quote = "", fill = T, comment.char = "!")

data_list = list(data3521, data10893, data14682, data14683, data37407, data63668)
GSE_list = list("GSE3521","GSE10893", "GSE14682", "GSE14683", "GSE37407","GSE63668")
GPL_list = list("GPL1390","GPL8128","GPL8128","GPL1390","GPL13703","GPL16686")

load("compare.Rdata")


# Get gene expression file
# Select samples from our data
for (i in c(1,2,5,6)){
  data = data_list[i][[1]]
  head(data)
  rownames(data) = data[,1]
  colnames(data) = substr(colnames(data),3,nchar(colnames(data))-1)
  data= data[,-1]
  validate = ifelse(compare$Dataset_id == as.character(GSE_list[i]),
                      as.character(compare$Sample_id),"")
  validate = subset(validate,validate != "")
  validate
  data = data[,c(validate)]
  head(data)
  data_list[i][[1]] = data
}

data3521 = data_list[1][[1]]
data10893 = data_list[2][[1]]
data37407 = data_list[5][[1]]
data63668 = data_list[6][[1]] 

validate = ifelse(compare$Dataset_id == as.character(GSE_list[3]),
                  as.character(compare$Sample_id),"")
validate = subset(validate,validate != "")
validate
data14682 = data14682[,c(validate)]
data14682 <- as.data.frame(data14682)

validate = ifelse(compare$Dataset_id == as.character(GSE_list[4]),
                  as.character(compare$Sample_id),"")
validate = subset(validate,validate != "")
validate
data14683 = data14683[,c(validate)]
data14683 = as.data.frame(data14683)

 

# Get platform information
library(GEOquery)
gpl1390 <- getGEO('GPL1390', destdir = ".")
gpl8128 <- getGEO('GPL8128', destdir = ".")
gpl13703 <- getGEO('GPL13703', destdir = ".")
gpl16686 <- getGEO('GPL16686',destdir = ".")

id1390 <- gpl1390@dataTable@table[,c("ID", "GENE_NAME")]
id8128 <- gpl8128@dataTable@table[,c("ID", "SYMBOL")]
id13703 <- gpl13703@dataTable@table[,c("Gene ID", "Gene Symbol")]
# Get gene name for GPL16686
anno16686 = read.table("refGene.txt", sep = "\t", header = F)
anno16686=anno16686[,c(2,13)]
names(anno16686)=c("acc","symbol")
anno16686=unique(anno16686)
rownames(anno16686)=anno16686$acc
symbol = as.character(anno16686[gpl16686@dataTable@table$GB_ACC, "symbol"])
symbol[is.na(symbol)] = ""
gpl16686@dataTable@table$symbol = symbol
id16686 <- gpl16686@dataTable@table[,c("ID", "symbol")]

# Filter the unique gene symbol
# Let gene symbol be the rownames
data3521 = data3521[rownames(data3521) %in% id1390$ID,]
id1390 = id1390[match(rownames(data3521), id1390$ID),]
tmp = by(data3521,
         id1390$GENE_NAME,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
data3521 = data3521[rownames(data3521) %in% probes, ]
rownames(data3521) = id1390[match(rownames(data3521),id1390$ID),2]

data10893 = data10893[rownames(data10893) %in% id1390$ID,]
id1390 = id1390[match(rownames(data10893), id1390$ID),]
tmp = by(data10893,
         id1390$GENE_NAME,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
data10893 = data10893[rownames(data10893) %in% probes, ]
rownames(data10893) = id1390[match(rownames(data10893),id1390$ID),2]

data14682 = data14682[rownames(data14682) %in% id8128$ID,]
id8128 = id8128[match(rownames(data14682), id8128$ID),]
tmp = by(data14682,
         id8128$SYMBOL,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
data14682 = data14682[rownames(data14682) %in% probes, ]
rownames(data14682) = id8128[match(rownames(data14682),id8128$ID),2]

data14683 = data14683[rownames(data14683) %in% id8128$ID,]
id8128 = id8128[match(rownames(data14683), id8128$ID),]
tmp = by(data14683,
         id8128$SYMBOL,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
data14683 = data14683[rownames(data14683) %in% probes, ]
rownames(data14683) = id8128[match(rownames(data14683),id8128$ID),2]

data37407 = data37407[rownames(data37407) %in% id13703$`Gene ID`,]
id13703 = id13703[match(rownames(data37407), id13703$`Gene ID`),]
tmp = by(data37407,
         id13703$`Gene Symbol`,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
data37407 = data37407[rownames(data37407) %in% probes, ]
rownames(data37407) = id13703[match(rownames(data37407),id13703$`Gene ID`),2]

data63668 = data63668[rownames(data63668) %in% id16686$ID,]
head(data63668)
id16686 = id16686[match(rownames(data63668), id16686$ID),]
tmp = by(data63668,
         id16686$symbol,
         function(x) rownames(x)[which.max(rowMeans(x))])
probes = as.character((tmp))
data63668 = data63668[rownames(data63668) %in% probes, ]
rownames(data63668) = id16686[match(rownames(data63668),id16686$ID),2]

save(data3521,data10893,data14682,data14683,data37407,data63668, file = "exp.Rdata")


# Get group information
rm(list = ls())
load("exp.Rdata")
load("compare.Rdata")
gl3521 = ifelse(compare[compare$Dataset_id == "GSE3521" &compare$Sample_id %in% colnames(data3521),"Sample_label"] == "Primary Tumor", 
                 "control", "treat")
gl10893 = ifelse(compare[compare$Dataset_id == "GSE10893" &compare$Sample_id %in% colnames(data10893),"Sample_label"] == "Primary Tumor", 
                 "control", "treat")
gl14682 = ifelse(compare[compare$Sample_id %in% colnames(data14682),"Sample_label"] == "Primary Tumor", 
                 "control", "treat")
gl14683 = ifelse(compare[compare$Sample_id %in% colnames(data14683),"Sample_label"] == "Primary Tumor", 
                 "control", "treat")
gl37407 = ifelse(compare[compare$Sample_id %in% colnames(data37407),"Sample_label"] == "Primary Tumor", 
                 "control", "treat")
gl63668 = ifelse(compare[compare$Sample_id %in% colnames(data63668),"Sample_label"] == "Primary Tumor", 
                 "control", "treat")
save(gl3521, gl10893, gl14682, gl14683, gl37407, gl63668, file = "group_list.Rdata")


# PCA
rm(list = ls())
load("exp.Rdata")
load("group_list.Rdata")
library(ggfortify)
library(ggpubr)
df3521 = as.data.frame(t(data3521))
df10893 = as.data.frame(t(data10893))
df14682 = as.data.frame(t(data14682))
df14683 = as.data.frame(t(data14683))
df37407 = as.data.frame(t(data37407))
df63668 = as.data.frame(t(data63668))

df3521$group = gl3521
df10893$group = gl10893
df14682$group = gl14682
df14683$group = gl14683
df37407$group = gl37407
df63668$group = gl63668

p1 <- autoplot(prcomp(df3521[,1:(ncol(df3521) - 1)]), data = df3521, colour = 'group', alpha = 0.6) + 
  labs(subtitle = "GSE3521") + 
  theme(axis.title = element_text(size = 9, face = "italic"), plot.subtitle = element_text(size = 9, hjust = 0.5))
p2 <- autoplot(prcomp(df10893[,1:(ncol(df10893) - 1)]), data = df10893, colour = 'group', alpha = 0.6) + 
  labs(subtitle = "GSE10893") + 
  theme(axis.title = element_text(size = 9, face = "italic"), plot.subtitle = element_text(size = 9, hjust = 0.5))
p3 <- autoplot(prcomp(df14682[,1:(ncol(df14682) - 1)]), data = df14682, colour = 'group', alpha = 0.6) + 
  labs(subtitle = "GSE14682") + 
  theme(axis.title = element_text(size = 9, face = "italic"), plot.subtitle = element_text(size = 9, hjust = 0.5))
p4 <- autoplot(prcomp(df14683[,1:(ncol(df14683) - 1)]), data = df14683, colour = 'group', alpha = 0.6) + 
  labs(subtitle = "GSE14683") + 
  theme(axis.title = element_text(size = 9, face = "italic"), plot.subtitle = element_text(size = 9, hjust = 0.5))
p5 <- autoplot(prcomp(df37407[,1:(ncol(df37407) - 1)]), data = df37407, colour = 'group', alpha = 0.6) + 
  labs(subtitle = "GSE37407") + 
  theme(axis.title = element_text(size = 9, face = "italic"), plot.subtitle = element_text(size = 9, hjust = 0.5))
p6 <- autoplot(prcomp(df63668[,1:(ncol(df63668) - 1)]), data = df63668, colour = 'group', alpha = 0.6) + 
  labs(subtitle = "GSE63668") + 
  theme(axis.title = element_text(size = 9, face = "italic"), plot.subtitle = element_text(size = 9, hjust = 0.5))


ggarrange(p1,p2,p3,p4,p5,p6, ncol = 3, nrow = 2, common.legend = T)


# Matrix preparation
rm(list = ls())
load("exp.Rdata")
load("group_list.Rdata")
options(stringsAsFactors = F)
library(limma)

design14682 <- model.matrix(~0+factor(gl14682))
colnames(design14682) = levels(factor(gl14682))
rownames(design14682) = colnames(data14682)
design14682

design14683 <- model.matrix(~0+factor(gl14683))
colnames(design14683) = levels(factor(gl14683))
rownames(design14683) = colnames(data14683)
design14683

design37407 <- model.matrix(~0+factor(gl37407))
colnames(design37407) = levels(factor(gl37407))
rownames(design37407) = colnames(data37407)
design37407

design63668 <- model.matrix(~0+factor(gl63668))
colnames(design63668) = levels(factor(gl63668))
rownames(design63668) = colnames(data63668)
design63668

design3521 <- model.matrix(~0+factor(gl3521))
colnames(design3521) = levels(factor(gl3521))
rownames(design3521) = colnames(data3521)
design3521


contrast.matrix <- makeContrasts(paste0(c("treat","control"),collapse = "-"),levels = design14682)
contrast.matrix


# Differentially expressed gene analysis with limma
fit14682 <- lmFit(data14682,design14682)
fit14682.2 <- contrasts.fit(fit14682,contrast.matrix)
fit14682.2 <- eBayes(fit14682.2)
tempOutput = topTable(fit14682.2, coef = 1, n=Inf)
nrDEG14682 = na.omit(tempOutput)
head(nrDEG14682)

fit14683 <- lmFit(data14683,design14683)
fit14683.2 <- contrasts.fit(fit14683,contrast.matrix)
fit14683.2 <- eBayes(fit14683.2)
tempOutput = topTable(fit14683.2, coef = 1, n=Inf)
nrDEG14683 = na.omit(tempOutput)
head(nrDEG14683)

fit37407 <- lmFit(data37407,design37407)
fit37407.2 <- contrasts.fit(fit37407,contrast.matrix)
fit37407.2 <- eBayes(fit37407.2)
tempOutput = topTable(fit37407.2, coef = 1, n=Inf)
nrDEG37407 = na.omit(tempOutput)
head(nrDEG37407)

fit63668 <- lmFit(data63668,design63668)
fit63668.2 <- contrasts.fit(fit63668,contrast.matrix)
fit63668.2 <- eBayes(fit63668.2)
tempOutput = topTable(fit63668.2, coef = 1, n=Inf)
nrDEG63668 = na.omit(tempOutput)
head(nrDEG63668)

# two channel analyzed with GEO2R (GSE3521 is identical with GSE10893 after selection)
# only perform analysis on GSE 3521
fit3521 <- lmFit(data3521, design3521) 
fit3521.2 <- contrasts.fit(fit3521, contrast.matrix)
fit3521.2 <- eBayes(fit3521.2)
tempOutput = topTable(fit3521.2, coef = 1, n=Inf)
nrDEG3521 = na.omit(tempOutput)
head(nrDEG3521)



save(nrDEG14682,nrDEG14683,nrDEG37407,nrDEG63668, nrDEG3521, file = "DEG.Rdata")

# Draw volcano plot
rm(list = ls())
load("DEG.Rdata")
library(ggplot2)

DEG_14682 = nrDEG14682
logFC_cf_14682 <- with(DEG_14682, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
DEG_14682$change = as.factor(ifelse(DEG_14682$P.Value < 0.05 & abs(DEG_14682$logFC) > logFC_cf_14682,
                                    ifelse(DEG_14682$logFC > logFC_cf_14682, "up", "down"), "not"))
title_14682 <- ("GSE14682: TREAT vs CONTROL")
sfg_14682 = DEG_14682[DEG_14682$change != "not", ]
g14682 = ggplot(data=DEG_14682, aes(x=logFC, y=-log10(P.Value), color=change, size = change)) +
  geom_point(alpha=1) +
  xlab("log2(fold change)") + ylab("-log10(Pvalue)") +
  ggtitle(title_14682) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
        plot.title = element_text(size=14,hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.box.background = element_rect(fill = NA, size = 0.8), 
        legend.position = c(0,0), legend.justification = c(0,0), legend.key.size = unit(3,"pt")) +
  scale_colour_manual(values = c('deepskyblue','black','red'), name = "p.value<0.05") +
  scale_size_manual(values = c(1.75, 0.5, 1.75)) +
  guides(size = FALSE)
g14682

DEG_14683 = nrDEG14683
logFC_cf_14683 <- with(DEG_14683, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
DEG_14683$change = as.factor(ifelse(DEG_14683$P.Value < 0.05 & abs(DEG_14683$logFC) > logFC_cf_14683,
                                    ifelse(DEG_14683$logFC > logFC_cf_14683, "up", "down"), "not"))
title_14683 <- ("GSE14683: TREAT vs CONTROL")
sfg_14683 = DEG_14683[DEG_14683$change != "not", ]
g14683 = ggplot(data=DEG_14683, aes(x=logFC, y=-log10(P.Value), color=change, size = change)) +
  geom_point(alpha=1) +
  xlab("log2(fold change)") + ylab("-log10(Pvalue)") +
  ggtitle(title_14683) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
        plot.title = element_text(size=14,hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.box.background = element_rect(fill = NA, size = 0.8), 
        legend.position = c(0,0), legend.justification = c(0,0), legend.key.size = unit(3,"pt")) +
  scale_colour_manual(values = c('black','red'), name = "p.value<0.05") +
  scale_size_manual(values = c(0.5, 1.75)) +
  guides(size = FALSE)
g14683



DEG_37407 = nrDEG37407
logFC_cf_37407 <- with(DEG_37407, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
DEG_37407$change = as.factor(ifelse(DEG_37407$P.Value < 0.05 & abs(DEG_37407$logFC) > logFC_cf_37407,
                                    ifelse(DEG_37407$logFC > logFC_cf_37407, "up", "down"), "not"))
title_37407 <- ("GSE37407: TREAT vs CONTROL")
sfg_37407 = DEG_37407[DEG_37407$change != "not", ]
g37407 = ggplot(data=DEG_37407, aes(x=logFC, y=-log10(P.Value), color=change, size = change)) +
  geom_point(alpha=1) +
  xlab("log2(fold change)") + ylab("-log10(Pvalue)") +
  ggtitle(title_37407) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
        plot.title = element_text(size=14,hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.box.background = element_rect(fill = NA, size = 0.8), 
        legend.position = c(0,0), legend.justification = c(0,0), legend.key.size = unit(3,"pt")) +
  scale_colour_manual(values = c('black'), name = "p.value<0.05") +
  scale_size_manual(values = c(0.5)) +
  guides(size = FALSE)
g37407

DEG_63668 = nrDEG63668
logFC_cf_63668 <- with(DEG_63668, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
DEG_63668$change = as.factor(ifelse(DEG_63668$P.Value < 0.05 & abs(DEG_63668$logFC) > logFC_cf_63668,
                                    ifelse(DEG_63668$logFC > logFC_cf_63668, "up", "down"), "not"))
title_63668 <- ("GSE63668: TREAT vs CONTROL")
sfg_63668 = DEG_63668[DEG_63668$change != "not", ]
g63668 = ggplot(data=DEG_63668, aes(x=logFC, y=-log10(P.Value), color=change, size = change)) +
  geom_point(alpha=1) +
  xlab("log2(fold change)") + ylab("-log10(Pvalue)") +
  ggtitle(title_63668) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
        plot.title = element_text(size=14,hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.box.background = element_rect(fill = NA, size = 0.8), 
        legend.position = c(0,0), legend.justification = c(0,0), legend.key.size = unit(3,"pt")) +
  scale_colour_manual(values = c('deepskyblue','black','red'), name = "p.value<0.05") +
  scale_size_manual(values = c(1.75,0.5,1.75)) +
  guides(size = FALSE)
g63668

DEG_3521 = nrDEG3521
logFC_cf_3521 <- with(DEG_3521, mean(abs(logFC),na.rm = TRUE) + 2*sd(abs(logFC),na.rm = TRUE))
DEG_3521$change = as.factor(ifelse(DEG_3521$P.Value < 0.05 & abs(DEG_3521$logFC) > logFC_cf_3521,
                                    ifelse(DEG_3521$logFC > logFC_cf_3521, "up", "down"), "not"))
title_3521 <- ("GSE3521: TREAT vs CONTROL")
sfg_3521 = DEG_3521[DEG_3521$change != "not", ]
g3521 = ggplot(data=DEG_3521, aes(x=logFC, y=-log10(P.Value), color=change, size = change)) +
  geom_point(alpha=1) +
  xlab("log2(fold change)") + ylab("-log10(P.value)") +
  ggtitle(title_3521) + 
  theme_bw() +
  theme(panel.grid = element_blank(), panel.border = element_rect(size = 1),
        plot.title = element_text(size=14,hjust = 0.5, face = "bold"),
        axis.title.x = element_text(size = 9), axis.title.y = element_text(size = 9),
        legend.box.background = element_rect(fill = NA, size = 0.8), 
        legend.position = c(0,0), legend.justification = c(0,0), legend.key.size = unit(3,"pt")) +
  scale_colour_manual(values = c('deepskyblue','black','red'), name = "p.value<0.05") +
  scale_size_manual(values = c(1.75, 0.5, 1.75)) +
  guides(size = FALSE)
g3521

write.csv(sfg_14682, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\deg analysis\\DEG\\sfg_14682.csv")
write.csv(sfg_14683, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\deg analysis\\DEG\\sfg_14683.csv")
write.csv(sfg_37407, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\deg analysis\\DEG\\sfg_37407.csv")
write.csv(sfg_63668, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\deg analysis\\DEG\\sfg_63668.csv")
# write.csv(sfg_3521, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\deg analysis\\DEG\\sfg_3521.csv")
write.csv(sfg_3521, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\deg analysis\\DEG\\sfg_3521_pvalue.csv")


save(sfg_14682,sfg_14683,sfg_37407,sfg_63668,sfg_3521, file = "significant gene.Rdata")

library(ggpubr)
g10893 = g3521 + ggtitle("GSE10893: TREAT vs CONTROL")
  
g = ggarrange(g3521, g10893, g14682, g14683, g37407, g63668, nrow = 2, ncol = 3)
g


# RobustRankAggreg package for integration
rm(list = ls())
load("significant gene.Rdata")
library(RobustRankAggreg)
set.seed(123)

# Throw away GSE37407 because it is miRNA-seq
UP_3521 = sfg_3521[order(sfg_3521$logFC,decreasing = T), ]
UP_14682 = sfg_14682[order(sfg_14682$logFC,decreasing = T), ]
UP_14683 = sfg_14683[order(sfg_14683$logFC, decreasing = T), ]
UP_63668 = sfg_63668[order(sfg_63668$logFC, decreasing = T), ]

geneUP3521 <- rownames(UP_3521)
geneUP14682 <- rownames(UP_14682)
geneUP14683 <- rownames(UP_14683)
geneUP63668 <- rownames(UP_63668)

glist_UP <- list(geneUP3521, geneUP14682, geneUP14683, geneUP63668)
freq_UP = as.data.frame(table(unlist(glist_UP)))
ag_UP = aggregateRanks(glist_UP)
ag_UP$Freq = freq_UP[match(ag_UP$Name,freq_UP$Var1),2]
head(ag_UP,15)
sfg_UP = ag_UP[ag_UP$Score <0.05, ]
sfg_UP
write.csv(sfg_UP, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\deg analysis\\DEG\\sfg_UP.csv")


DOWN_3521 = sfg_3521[order(sfg_3521$logFC), ]
DOWN_14682 = sfg_14682[order(sfg_14682$logFC), ]
DOWN_14683 = sfg_14683[order(sfg_14683$logFC), ]
DOWN_63668 = sfg_63668[order(sfg_63668$logFC), ]

geneDOWN3521 <- rownames(DOWN_3521)
geneDOWN14682 <- rownames(DOWN_14682)
geneDOWN14683 <- rownames(DOWN_14683)
geneDOWN63668 <- rownames(DOWN_63668)

glist_DOWN <- list(geneDOWN3521, geneDOWN14682, geneDOWN14683, geneDOWN63668)
freq_DOWN = as.data.frame(table(unlist(glist_DOWN)))
ag_DOWN = aggregateRanks(glist_DOWN)
ag_DOWN$Freq = freq_DOWN[match(ag_DOWN$Name,freq_DOWN$Var1),2]
head(ag_DOWN,15)
sfg_DOWN = ag_DOWN[ag_DOWN$Score <0.05, ]
sfg_DOWN
write.csv(sfg_DOWN, file = "D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\deg analysis\\DEG\\sfg_DOWN.csv")



















