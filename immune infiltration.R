setwd("D:\\STUDY\\szbl\\metastatic cell genomic\\code\\Primary_CNS")
data <- read.csv("D:\\STUDY\\szbl\\metastatic cell genomic\\results\\Primary_CNS\\immune filtration\\estimation_matrix.csv", 
                  header = T)
library(ggplot2)
library(cowplot)
library(patchwork)
TIMER <- data[1:6, ]
rownames(TIMER) = TIMER$cell_type
TIMER = TIMER[,c(-1)]
load("group_list.Rdata")
gl <- c(gl14682, gl14683)

TIMER = as.data.frame(t(TIMER))
TIMER$group = gl
TIMER = TIMER[order(TIMER$group), ]
TIMER = TIMER[, -c(ncol(TIMER))]
TIMER = as.data.frame(t(TIMER))

gl = sort(gl)

# Plot pie
pie_plot <- function(pc){
  sample = as.data.frame(TIMER[,pc])
  rownames(sample) = substr(rownames(TIMER),1, nchar(rownames(TIMER))-6) 
  colnames(sample)[1] = colnames(TIMER)[pc]
  subtitle = paste(colnames(sample)[1], "-", gl[[pc]])
  p = ggplot(sample, aes(x = "", y = sample[,1], fill = rownames(sample))) +
    geom_bar(stat = "identity", width = 1) +
    coord_polar(theta = "y") +
    labs(x = "", y ="", title = "", subtitle = subtitle) +
    scale_fill_brewer(palette = "Set3") +
    theme(axis.ticks = element_blank(), axis.text.x = (element_blank()),
          legend.title = element_blank(),legend.position = "none",
          panel.grid = element_blank(), panel.border = element_blank(),
          plot.subtitle = element_text(hjust = 0.5, size = 8)
    )
}

pie_list <- lapply(1:80, pie_plot) # R might be collapsed


pp_patchwork <- wrap_plots(pie_list, byrow = T, nrow = 2) + plot_layout(guides = 'collect')

pp_patchwork


# Heatmap
options(stringsAsFactors = F)
library(pheatmap)
choose_gene = rownames(TIMER)
choose_matrix = TIMER[choose_gene,]
rownames(choose_matrix) = substr(rownames(TIMER),1, nchar(rownames(TIMER))-6)
matrix_scale = t(scale(t(choose_matrix)))
annotation_col = data.frame(sampletype = gl)
row.names(annotation_col) <- colnames(choose_matrix)
pheatmap(matrix_scale, scale = "row", annotation_col = annotation_col,
         cluster_rows = F, cluster_cols = T,show_rownames = T,show_colnames = F,border_color = NA)


