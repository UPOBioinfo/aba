###########################################
# Eugenio Mangas / Antonio J. Pérez, 2019 #
# Enrichment of core and accesory genome  #
###########################################

library(ggpubr)
library(topGO)
library(dplyr)
library(ggplot2)
library(cowplot)

palette <- c("#F52A2A", "#D561EA", "#61B0EA", "green", "#E89B57", "#E4EA61", "white")

# files
files <- c("g99-100.gn", "g20-99.gn", "g1-20.gn", "g0-1.gn")
file_bg <- "../pangenome/pangenome_annot2.tsv"
Nodes <- 21 # number of processes to show
Ontology <- "GO.P.ID" #PFC (BP MF CC)

#Create temp file
data <- read.csv(file_bg, sep = "\t", header = TRUE, row.names = NULL)[,(c('GENENAME', Ontology))] #other id: X.ID
data$GO.P.ID <- as.character(gsub(';', ', ', data$GO.P.ID))
file_temp <- paste0(file_bg,"2")
write.table(data, file = file_temp, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

###############################
# Iterate the different files #
###############################
n = 0
figure <- list()
mylabels <- list()
f <- list()

for (file_gene in files) {

n = n + 1
file_gene <- files[[n]]
# Get gene IDs for the enrichment
genes <- read.csv(file_gene, header=F)$V1

# Get background annotation
GOesByID <- readMappings(file = file_temp)
bg_genes <- names(GOesByID)

compared_genes <- factor(as.integer(bg_genes %in% genes))
names(compared_genes) <- bg_genes

# Create topGO object
GOdata <- new("topGOdata", ontology = "BP", allGenes = compared_genes,
              annot = annFUN.gene2GO, gene2GO = GOesByID)
asd<-unlist(Term(GOTERM))

# Run Fisher test
resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")

# Create and print table with enrichment result
allRes <- GenTable(GOdata, classicFisher = resultFisher, topNodes = Nodes)

# Graphic
#########
layout(t(1:2), widths=c(8,3))
par(mar=c(4, .5, .7, .7), oma=c(3, 15, 3, 4), las=1)

pvalue <- as.numeric(gsub("<", "", allRes$classicFisher)) # remove '<' symbols
allRes$classicFisher <- pvalue
max_value <- as.integer(max(-log(pvalue)))+1
pv_range <- exp(-seq(max_value, 0, -1))
allRes = mutate(allRes, plot_id = paste(GO.ID, Term, sep = " - "))

mylabels[[n]] <- paste (allRes$GO.ID, "-",  asd[allRes$GO.ID])
mybreaks = c(1,1.0e-5,1.0e-10,1.0e-15,1.0e-20,1.0e-25,1.0e-30,1.0e-40,0)

figure[[n]] <- ggplot(data=allRes, aes(x=reorder(plot_id, Significant), y=Significant) ) +
  geom_bar(stat="identity", color="black", aes(fill=as.numeric(log(classicFisher))), size = 0.3)+
  geom_text(aes(label=mylabels[[n]]), position=position_fill(vjust=0), hjust=0, fontface="bold", size = 5)+
  coord_flip() +
  theme(panel.background = element_blank(), panel.grid.major.x = element_line(colour = "darkgrey", size=0.75),
        panel.grid.minor.x = element_line(colour = "grey",size=0.75), axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        axis.ticks.x =element_blank(), axis.line.y=element_blank(), 
        legend.position = "none") +
  ylab("number of genes") + 
  guides(fill = guide_colourbar(barheight = 25, reverse=T)) +
  scale_fill_gradientn(name = "p-value", colours = palette, limits = log(c(1e-31,1)), breaks = log(mybreaks), 
                   guide = guide_colourbar(reverse = TRUE), labels=mybreaks)
f[[n]] <- ggarrange(figure[[n]])

}

# Final figure
empty <-
  ggplot(data=allRes, aes(x=reorder(plot_id, Significant), y=Significant) ) +
  geom_bar(stat="identity", color="black", aes(fill=as.numeric(log(classicFisher))), size = 0.3)+
  theme(panel.background = element_blank(), panel.grid.major.x = element_line(colour = "darkgrey",size=0.75),
        panel.grid.minor.x = element_line(colour = "grey",size=0.5), axis.title.y=element_blank(), 
        axis.title.x=element_blank(), axis.text.y=element_blank(), axis.ticks.y=element_blank(), 
        axis.ticks.x =element_blank(), axis.line.y=element_blank(), axis.line.x=element_blank())+
  ylab("number of genes") + 
  guides(fill = guide_colourbar(barheight = 25, barwidth = 1.6, reverse=T)) +
  scale_fill_gradientn(name = "p-value", colours = palette, limits = log(c(1e-31,1)), breaks = log(mybreaks), 
                       guide = guide_colourbar(reverse = TRUE), labels=mybreaks)
my_legend <- get_legend(empty)
legend <- as_ggplot(my_legend)

pdf("../figures/Fig3.pdf", width=16, height=8, paper='special', onefile=FALSE)
figab <- ggarrange(f[[1]], f[[2]], labels = c("a)", "b)"), font.label = list(size = 16, color = "black", face = "bold"))
figcd <- ggarrange(f[[3]], f[[4]], labels = c("c)", "d)"), font.label = list(size = 16, color = "black", face = "bold"))
figabcd <- ggarrange(figab, figcd, nrow = 2)
ggarrange(figabcd, legend, font.label = list(size = 14, color = "black", face = "bold", family = NULL),
          common.legend = TRUE, legend = "right", nrow = 1, widths = c(1, 0.1))
dev.off()

