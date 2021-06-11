# Check both number of genes and shared genes by strain
# Change working directory, pan_clusters2_table.matrix and metadata.tsv
#  The latter file should have 5 columns, though the 4 latter are empty: Number, isolation, whyDeleted, COLLECTION_DATE, Proteins_Prokka
# AJPerez, 2019 (updated June 2021)
library(magrittr)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(ggrepel)
library(ggpubr)

setwd("./")

# Palette
c10 <- c("red","green","darkorange","darkturquoise","darkred","grey",
         "magenta","darkgreen","blue","yellow3")
s4 <- c(1,2,4,3)

# Gather data
data <- read.table("pan_clusters2_table.matrix", comment.char="#", sep="\t")
rnames <- data[,1]                          # assign labels in column 1 to "rnames"
pgenes <- data.matrix(data[,2:ncol(data)])  # transform columns 2-> into a matrix
rownames(pgenes) <- rnames                  # assign row names
colnames(pgenes) <- rnames

# Metadata
mdata <- read.csv("metadata.tsv", sep="\t")
isolation <- mdata$"isolation"
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
  geom_point(aes(), size=1, stroke=2.5) +
  scale_colour_manual(values = c10) +
  scale_shape_manual(values = s4) +
  xlab("Number of genes") +
  ylab("Average number of shared genes") +
  theme(plot.title = element_text(hjust = -0.05), text = element_text(size=16), legend.text = element_text(size = 13),
        legend.position = "top", legend.box = "horizontal", legend.title = element_blank(), 
        legend.spacing.x = unit(0.1, 'cm')) +
  guides(color = guide_legend(nrow = 2, keyheight = 1.3), shape =  guide_legend(nrow = 2, keyheight = 1.3)) +
  scale_x_continuous(breaks = seq(0, 5000, by = 500), limits =c(0, 5000))
p

length2 <- data.frame(n=x, total=mdata$Proteins_Prokka)
length2 <- stack(length2)

redline1 <- mean(length2$values) - sd(length2$values) * 2.5
h1 <- ggplot(data=length2, aes(values, fill=ind)) + 
  geom_histogram(binwidth = 1, color="grey20", aes(fill=ind), position = "dodge") +
  xlab("Number of genes") +
  ylab("Number of strains") +
  scale_x_continuous(limits =c(3250, 4000)) +
  theme(plot.title = element_text(hjust = -0.1)) +
  geom_density(aes(y=1*..count..), size = 1, alpha = .5) +
  geom_vline(xintercept = redline1, color="red", linetype="dashed", size=1.25) +
  theme_grey() +
  theme(legend.title = element_blank(), legend.position=c(0.82, 0.87), text = element_text(size=16))
h1

bluepoint <- 3136
d2 <- data_length %>% filter(shared < bluepoint)
redline2 <- mean(d2$shared) - sd(d2$shared) * 2.5
h2 <- ggplot(data=data_length, aes(shared, label=id)) + 
  geom_histogram(binwidth = 1, color="grey40", fill="grey") +
  xlab("Average number of shared genes") +
  ylab("Number of strains") +
  scale_x_continuous(breaks = seq(0, 4000, by = 50), limits =c(1500, 1700)) +
  geom_density(aes(y=1*..count..), size = 1.2) +
  geom_vline(xintercept = redline2, color="red", linetype="dashed", size=1.25) +
  geom_vline(xintercept = bluepoint, color="blue", linetype="dashed", size=1.25) +
  theme_grey() +
  theme(plot.title = element_text(hjust = -0.1), text = element_text(size=16))
h2

diptest::dip.test(length2$values)
diptest::dip.test(data_length$shared)

ab <- ggarrange(h1, h2, labels = c("a)", "b)"), font.label = list(size = 16, color = "black", face = "bold"), nrow = 2)

pdf("ngenes.pdf", width=16, height=8, paper='special')
ggarrange(ab, p, labels = c("", "c)"), font.label = list(size = 16, color = "black", face = "bold"),  widths = c(1/4, 3/4))
dev.off()

# Groups/Peaks
#g1 <- which(data_length$shared < 3136)
#write.csv(g1, file = "g1.tsv")
#g2 <- which(data_length$shared >= 3136)
#write.csv(g2, file = "g2.tsv")

