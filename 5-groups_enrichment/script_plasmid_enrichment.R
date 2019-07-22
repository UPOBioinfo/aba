###########################################
# Antonio J. Pérez, 2019                  #
# Enrichment of gene groups               #
###########################################

library(ggplot2)

data <- read.table("group12_enriched.tsv", sep = "\t")

min <- min(data$V5)
p1 <- data.frame (xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=0)
p2 <- data.frame (xmin=-Inf, xmax=Inf, ymin=0, ymax=Inf)
p <- ggplot() +
  geom_bar(data=data, aes(x = reorder(V1, V4), y = -V2, fill = V3), stat = "identity") +
  geom_rect(data=p2, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill="blue", alpha=0.1, inherit.aes = FALSE) +
  geom_bar(data= data, aes(x = reorder(V1, V4), y = V4, fill = V5), stat = "identity") +
  geom_rect(data=p1, aes(xmin=xmin, xmax=xmax, ymin=ymin, ymax=ymax), 
            fill="green", alpha=0.1, inherit.aes = FALSE) +
  geom_bar(data=data, aes(x = reorder(V1, V4), y = -V2, fill = V3), stat = "identity") +
  scale_fill_gradient(low="darkred", high="#F1C9C9",
                      trans="log10",name = "p-value",
                      breaks=c(1.0e-18, 1.0e-16, 1.0e-14, 1e-12, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2, 1)) +
  coord_flip() + 
  theme(text=element_text(size=16), axis.text = element_text(color = "black")) +
  ylab("Number of genes") +
  xlab("Functional annotations") +
  guides(fill = guide_colourbar(barheight = 25, reverse=T)) +
  annotate("text", x=0.83, y=58, label= "group 2", color="darkblue", size=5) +
  annotate("text", x=0.83, y=-26, label= "group 1", color="darkgreen", size=5) +
  geom_hline(yintercept = 0, color = "white") +
  scale_y_continuous(breaks = c(-40, -30, -20, -10, 0, 10, 20, 30, 40, 50, 60, 70),
                     labels = c(40, 30, 20, 10, 0, 10, 20, 30, 40, 50, 60, 70))
p

pdf("../figures/Fig4.pdf", width=10, height=6, paper='special')
p
dev.off()
