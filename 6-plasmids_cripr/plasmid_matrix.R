########################
# Eugenio Mangas, 2019 #
# Matrix of plasmids   #
########################

library(pheatmap)
library(viridisLite)
library(viridis)
library(dendextend)

mtrx<-read.csv("ab_plasm_matrix_2112_84_48without0_v3.tsv", sep = "\t", header = TRUE, row.names = "X")

mtrx[is.na(mtrx)] = 0
mtrx<-as.matrix(mtrx)
mtdt <- read.csv("groups_strains2_v3.tsv", check.names = FALSE, sep = "\t", header = TRUE, row.names = rownames(mtrx))
mtrx_cor <- cor(mtrx)


annotation_colors = list("Group" = c("Group 1 (Cas genes)"="purple","Group 1 (no-Cas genes)"="magenta","Group 2"="chartreuse2"))

#pheatmap(mtrx_cor, viridis(10),cluster_cols = TRUE, cluster_rows = TRUE)

pdf("../figures/SupplFigure4.pdf", width=10, height=6, paper='special')
pheatmap(mtrx, color = viridis(2), show_rownames=F, treeheight_row = 0, legend = FALSE, angle_col = 90, fontsize = 13,
         cluster_cols = TRUE, cluster_rows = FALSE, annotation_row=mtdt, annotation_colors = annotation_colors)
dev.off()

