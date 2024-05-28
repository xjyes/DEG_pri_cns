setwd("D:\\STUDY\\szbl\\metastatic cell genomic\\code")

# Load enriched gene alternation
data = read.csv("D:\\STUDY\\szbl\\metastatic cell genomic\\dataset\\Primary_CNS\\Primary_CNS_combined.csv")
all = read.csv("D:\\STUDY\\szbl\\metastatic cell genomic\\dataset\\Primary_CNS\\Primary_CNS_all.csv")

library(dplyr)
library(tidyr)
data = separate(data = data, col = 3 ,into = c("Primary#","Primary%"),sep = "[(]" )
data$`Primary%` = substr(data$`Primary%`,1,nchar(data$`Primary%`)-2)
data = separate(data = data, col = 5 ,into = c("CNS#","CNS%"),sep = "[(]" )
data$`CNS%` = substr(data$`CNS%`,1,nchar(data$`CNS%`)-2)
colnames(data)[1] = "Gene"

all = separate(data = all, col = 3 ,into = c("Primary#","Primary%"),sep = "[(]" )
all$`Primary%` = substr(all$`Primary%`,1,nchar(all$`Primary%`)-2)
all = separate(data = all, col = 5 ,into = c("CNS#","CNS%"),sep = "[(]" )
all$`CNS%` = substr(all$`CNS%`,1,nchar(all$`CNS%`)-2)
write.csv(all,file = "all.csv", row.names = F)

amplification <- subset(data,data$Alternation.Type == "Amplification")
deletion <- subset(data,data$Alternation.Type == "Deletion")
mutation <- subset(data,data$Alternation.Type == "Mutation")

p_amplification = data.frame(Percentage = c(amplification$`Primary%`,amplification$`CNS%`),
                             Sample.Type = c(rep("Primary",nrow(amplification)),rep("CNS",nrow(amplification))),
                             Gene.Symbol= c(rep(amplification$Gene,2)))
p_deletion = data.frame(Percentage = c(deletion$`Primary%`,deletion$`CNS%`),
                             Sample.Type = c(rep("Primary",nrow(deletion)),rep("CNS",nrow(deletion))),
                             Gene.Symbol= c(rep(deletion$Gene,2)))
p_mutation = data.frame(Percentage = c(mutation$`Primary%`,mutation$`CNS%`),
                             Sample.Type = c(rep("Primary",nrow(mutation)),rep("CNS",nrow(mutation))),
                             Gene.Symbol= c(rep(mutation$Gene,2)))

# Draw boxplot by the alteration type

library(ggplot2)
library(ggpubr)
b_a <- ggplot(p_amplification, aes(x = Sample.Type, y = as.numeric(Percentage))) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.5, position=position_dodge(0.6)) +
  geom_boxplot(show.legend = TRUE, width = 0.5, fill = "#A9A9A9") +
  geom_point(aes(color = Gene.Symbol),size = 2) +
  geom_line(aes(group = Gene.Symbol, color = Gene.Symbol), lwd = 0.5) +
  geom_text(aes(x = Sample.Type, y = as.numeric(Percentage)), label = p_amplification$Percentage, size = 3) +
  theme(panel.grid = element_line(linetype = "dotted",color = "grey"), axis.line = element_line(colour = "black", size = 1),
        panel.background = element_blank(),plot.title = element_text(size = 15,hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = 10, color = "black"),axis.title = element_text(size = 12, color = "black"))+
  labs(x = '', y = '% of samples that have an alteration', title = 'Amplification in Breast Cancer', subtitle = 'q-value<0.05')
b_a


b_d <- ggplot(p_deletion, aes(x = Sample.Type, y = as.numeric(Percentage))) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.5, position=position_dodge(0.6)) +
  geom_boxplot(show.legend = TRUE, width = 0.5, fill = "#A9A9A9") +
  geom_point(aes(color = Gene.Symbol),size = 2) +
  geom_line(aes(group = Gene.Symbol, color = Gene.Symbol), lwd = 0.5) +
  geom_text(aes(x = Sample.Type, y = as.numeric(Percentage)), label = p_deletion$Percentage, size = 3) +
  theme(panel.grid = element_line(linetype = "dotted",color = "grey"), axis.line = element_line(colour = "black", size = 1),
        panel.background = element_blank(),plot.title = element_text(size = 15,hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = 10, color = "black"),axis.title = element_text(size = 12, color = "black"))+
  labs(x = '', y = '% of samples that have an alteration', title = 'Deletion in Breast Cancer', subtitle = 'q-value<0.05')
b_d

b_m <- ggplot(p_mutation, aes(x = Sample.Type, y = as.numeric(Percentage))) +
  stat_boxplot(geom = "errorbar", width = 0.3, size = 0.5, position=position_dodge(0.6)) +
  geom_boxplot(show.legend = TRUE, width = 0.5, fill = "#A9A9A9") +
  geom_point(aes(color = Gene.Symbol),size = 2) +
  geom_line(aes(group = Gene.Symbol, color = Gene.Symbol), lwd = 0.5) +
  geom_text(aes(x = Sample.Type, y = as.numeric(Percentage)), label = p_mutation$Percentage, size = 3) +
  theme(panel.grid = element_line(linetype = "dotted",color = "grey"), axis.line = element_line(colour = "black", size = 1),
        panel.background = element_blank(),plot.title = element_text(size = 15,hjust = 0.5),
        plot.subtitle = element_text(size = 12, hjust = 0.5),
        axis.text = element_text(size = 10, color = "black"),axis.title = element_text(size = 12, color = "black"))+
  labs(x = '', y = '% of samples that have an alteration', title = 'Mutation in Breast Cancer', subtitle = 'q-value<0.05')
b_m

b_combined <- ggarrange(plotlist = list(b_a,b_d,b_m),
          ncol = 3, nrow = 1)
b_combined
ggsave(filename = "Gene alternation in Breast Cancer.png",
       width = 14, height = 8, unit = "in")

# Draw bar plot for all significant gene alteration
library(reshape2)
library(grid)
all <- melt(all, measure.vars = c("Primary%","CNS%"),variable.name = "Groups",value.name = "Percentage")
all$Gene<- factor(all$Gene,levels=c(all$Gene[(nrow(all)/2):1]),order = TRUE)


text_CNS <- textGrob("CNS -->", gp=gpar(fontsize=8, fontface="italic"))
text_Primary <- textGrob("<-- Primary", gp=gpar(fontsize=8, fontface="italic"))

p <- ggplot(all, aes(x=Gene, y=ifelse(Groups=="CNS%", as.numeric(Percentage),as.numeric(Percentage)*-1))) +
  annotation_custom(text_CNS,xmin = -1.5, xmax = -1.5, ymin = 8,ymax = 8) +
  annotation_custom(text_Primary,xmin = -1.5,xmax = -1.5, ymin = -10,ymax = -10) +
  scale_y_continuous(expand=expansion(add = c(0.1, 0.1)),
                     limits = c(-40, 80),
                     breaks = c(-25, 0,25,50,75),
                     label = c("25", "0","25","50","75")) +
  coord_flip(clip = "off") +
  theme(text=element_text(family = "sans",colour ="black",size = 12),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(vjust = -2),
        axis.line.y=element_blank(),
        axis.line = element_line(size = 1,colour = "black"),
        axis.ticks = element_line(size =1,colour = "black"),
        axis.ticks.length = unit(1.5,units = "mm"),
        panel.grid = element_line(linetype = "dotted",color = "grey"),
        panel.background = element_blank(),
        plot.margin=unit(x=c(0.2,0.1,0.5,0.1),
                          units="inches")) +
  geom_bar(aes(fill = q.Value), color = "white", width = 0.98,size = 1,stat = "identity") +
  geom_text(aes(label=ifelse(Groups == "CNS%",Percentage,"")),position = position_dodge2(width = 0.9,reverse = T),
              vjust = 0.8, hjust = 0, size = 2.5)+
  geom_text(aes(label=ifelse(Groups == "Primary%",Percentage,"")),position = position_dodge2(width = 0.9,reverse = T),
            vjust = 0, hjust = 1,  size = 2.5)+
  scale_fill_gradient(low = "green",high = "red") +
  labs(fill="q.Value",x="",y="% of samples altered",title="Gene alternation in Breast Cancer")


p
ggsave(filename = "Gene alternation in Breast Cancer_summary.png",
       width = 6, height = 6, unit = "in")









