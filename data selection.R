setwd("D:\\STUDY\\szbl\\metastatic cell genomic\\code")

data = read.csv("D:\\STUDY\\szbl\\metastatic cell genomic\\dataset\\dataset_information.csv")
head(data)

library(tidyr)

data = separate(data, col = Standard, into = c("comparison","site"), sep = "]",remove = F)
data["comparison"] = substr(data$comparison,2,nchar(data$comparison))

sort(table(data$comparison),decreasing = T)

brain_mets <- subset(data,data$Metastasis_site == "brain")

sort(table(brain_mets$Dataset_id),decreasing = T)
sort(table(brain_mets$comparison),decreasing = T)
sort(table(brain_mets$site),decreasing = T)
sort(table(brain_mets$Standard),decreasing = T)

compare <- data[data$comparison == "primary tumor vs. metastasis tumor"&
                  data$Metastasis_site == "brain",]

table(compare[compare$Dataset_id == "GSE14682","Platform_id"])
save(compare, file = "compare.Rdata")
