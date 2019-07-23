library(magrittr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)

setwd("/home/ajperez/Dropbox/aba/pangenome")

# Palette
c10 <- c("red","green","darkorange","darkturquoise","darkred","grey",
         "magenta","darkgreen","blue","yellow3")
s4 <- c(1,2,4,3)

# Gather data
data <- read.table("matrix_gn.tsv", comment.char="#", sep="\t")
rnames <- data[,1]                          # assign labels in column 1 to "rnames"
pgenes <- data.matrix(data[,2:ncol(data)])  # transform columns 2-> into a matrix
rownames(pgenes) <- rnames                  # assign row names
colnames(pgenes) <- rnames

# Metadata
mdata <- read.csv("metadata_file_2112.tsv", sep="\t")
isolation <- mdata$Isol.Source.4
isolation2 <- as.numeric(as.factor(isolation))
legend.cols <- as.numeric(as.factor(levels(isolation)))
whydeleted <- mdata$whyDeleted
mdata <- mdata %>% filter(Number == rnames) # only strains in input


# Gather mean length by every strain pair
ngenes <- vector()
means <- matrix("", length(rownames(pgenes)), length(colnames(pgenes)))
for (i in 1:length(rownames(pgenes))) {
  ngenes[i] <- pgenes[i,i]
}

# Means by strain
rowmeans <- vector()
for (i in 1:length(rownames(pgenes))) {
  rowmeans[i] <- mean(pgenes[i,-i]) 
}

# Convert to Vectors
x <- ngenes   # means
y <- rowmeans # means

# ggplot plot
data_length <- data.frame(id=mdata$N, ngenes=x, shared=y, category=isolation, 
                          color=isolation2, whydeleted=whydeleted, year=mdata$COLLECTION_DATE, 
                          proteins=mdata$Proteins_Prokka)

p <- ggplot(data=data_length, aes(ngenes, shared, label=id))+ 
  theme_gray() +
  stat_smooth(method = "lm", lty = 2, col = "black") +
  geom_point(aes(col=category, shape=whydeleted), size=5, stroke=2.5) +
  scale_colour_manual(values = c10) +
  scale_shape_manual(values = s4) +
  xlab("Number of genes") +
  ylab("Average number of shared genes") +
  theme(plot.title = element_text(hjust = -0.05), text = element_text(size=16), legend.text = element_text(size = 13),
        legend.position = "top", legend.box = "horizontal", legend.title = element_blank()) +
  guides(color = guide_legend(nrow = 2, keyheight = 1.3), shape =  guide_legend(nrow = 2, keyheight = 1.3)) +
  scale_x_continuous(breaks = seq(1500, 4500, by = 500), limits =c(1670, 4450)) +
  scale_y_continuous(breaks = seq(0, 3500, by = 1000), limits =c(0, 3500))
p

length2 <- data.frame(shared=x, total=mdata$Proteins_Prokka)
length2 <- stack(length2)

h1 <- ggplot(data=length2, aes(values, fill=ind)) + 
  #geom_histogram(binwidth = 20, color="grey20", aes(fill=ind), position = "dodge") +
  xlab("Number of genes") +
  ylab("Number of strains") +
  scale_x_continuous(limits =c(2500, 4500)) +
  theme(plot.title = element_text(hjust = -0.1)) +
  geom_density(aes(y=20*..count..), size = 1, alpha = .5) +
#  geom_vline(xintercept = 3202, color="red", linetype="dashed", size=1.25) +
  theme_grey() +
  theme(legend.title = element_blank(), legend.position=c(0.82, 0.87), text = element_text(size=16))
h1

h2 <- ggplot(data=data_length, aes(shared, label=id)) + 
  geom_histogram(binwidth = 5, color="grey40", fill="grey") +
  xlab("Average number of shared genes") +
  ylab("Number of strains") +
  scale_x_continuous(breaks = seq(0, 4000, by = 500), limits =c(2000, 3600)) +
  theme(plot.title = element_text(hjust = -0.1), text = element_text(size=16)) + 
  geom_density(aes(y=5*..count..), size = 1.2) +
#  geom_vline(xintercept = 2588, color="red", linetype="dashed", size=1.25) +
#  geom_vline(xintercept = 3136, color="blue", linetype="dashed", size=1.25) +
  theme_grey() +
  theme(legend.title = element_blank(), legend.position=c(0.82, 0.87), text = element_text(size=16))
h2

ab <- ggarrange(h1, h2, labels = c("a)", "b)"), font.label = list(size = 16, color = "black", face = "bold"), nrow = 2)

pdf("../figures/SupplFigure1.pdf", width=16, height=8, paper='special')
ggarrange(ab, p, labels = c("", "c)"), font.label = list(size = 16, color = "black", face = "bold"),  widths = c(1/4, 3/4))
dev.off()

# Groups/Peaks
#g1 <- which(data_length$shared < 3136)
#write.csv(g1, file = "g1.tsv")
#g2 <- which(data_length$shared >= 3136)
#write.csv(g2, file = "g2.tsv")

