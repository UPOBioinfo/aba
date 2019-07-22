#########################
# Antonio J. Pérez 2019 #
# Plasmid distribution  #
#########################

library(ggplot2)
library(tidyverse)
library(ggrepel)
library(plotly)
library(reshape2)
library(ggnewscale)
library(ggbeeswarm)
library(ggpubr)

data2 <- read.csv(file="plasmids_crispr_0_100_v3.tsv", header=TRUE, sep="\t")
data <- read.csv(file="plasmids_group1_0_100_v3.tsv", header=TRUE, sep="\t")
data3 <- read.csv(file="plasmids_crispr_0_100_cas9_v3_only.tsv", header=TRUE, sep="\t")

data$total_plasmids <- as.factor(data$total_plasmids)

point_size <- 3
color2_2 = "#0000FF"
color1_2 = rgb(1,0,0,0.5)
color1_1 = "#E3E3E3"
color2_1 = "#E3E3E3"

n_fun <- function(x){
  return(data.frame(y = 2.25, size = 4,
                    label = paste0("n=", length(x))))
}

mm <- 'median'
xmm <- 1.25

#data2 %>% filter(set == "group 1") %>% summarise(mean = mean(nplasmids))
#data2 %>% filter(set == "group 2") %>% summarise(mean = mean(nplasmids))
#data3 %>% filter (prokka_arrays > 0 & set == "group 1")  %>% summarise(mean = mean(prokka_arrays), stdev = sd(prokka_arrays))
#data3 %>% filter (set == "group 1")  %>% summarise(mean = mean(prokka_arrays), stdev = sd(prokka_arrays))
#data2 %>% filter (set == "group 1")  %>% summarise(mean = mean(nplasmids), stdev = sd(nplasmids))

# error bars
#stat_summary(fun.y=mean, geom="point", position="dodge", color = "grey35") +
#stat_summary(fun.data=mean_sdl, geom="errorbar", position="dodge", width = 0.1, color = "grey35")
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  if (ymin < 0) {ymin=0}
  return(c(y=m,ymin=ymin,ymax=ymax))
}
meanbar <- stat_summary(fun.y=mean, geom="point", color="grey35")
errorbar <- stat_summary(fun.data=data_summary, geom="errorbar", color="grey35", width=0.1)

p3 <-
  ggplot(data2, aes(set, nplasmids)) + 
  geom_quasirandom(aes(color=crispr), alpha=1, size=point_size) +
  scale_color_gradient(name = "CRISPR arrays (plasmids)", low = color1_1, high = "red", breaks = seq(0, max(data2$crispr), 1), 
                       limits = c(0,max(data2$crispr))) +
  stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = mm, geom = 'crossbar', color = 'black', 
               linetype = "dashed") +
  theme_classic()  +
  xlab("all strains") +
  ylab("Number of plasmids") +
  scale_y_continuous(breaks = seq(0, max(data$nplasmids), by = 1), limits = c(0, max(data2$nplasmids))) +
  scale_x_discrete(label=c("Cas strains", "no-Cas strains")) +
  theme(legend.key.height = unit(0.4, "cm"), legend.key.width = unit(1.1, "cm"),
        text = element_text(size=12), legend.position = "none") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = xmm, vjust = -5,
               aes(group=set), position = position_dodge(0.6)) +
  errorbar + meanbar
p3

p4 <-
  data2 %>%
  filter (crispr > 0) %>%
  ggplot(aes(set, nplasmids)) + 
  geom_quasirandom(aes(color=crispr), alpha=1, size=point_size) +
  scale_color_gradient(name = "CRISPR arrays (plasmids)", low = color1_1, high = "red", breaks = seq(0, max(data2$crispr), 1),
                       limits = c(0,max(data2$crispr)), guide = guide_colourbar(barwidth=1, barheight=8)) +
  stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = mm, geom = 'crossbar', color = 'black', linetype = "dashed") +
  theme_classic()  +
  xlab("CRISPR arrays in plasmids > 0") +
  ylab("Number of plasmids") +
  labs(y=NULL) +
  scale_x_discrete(label=c("Cas strains", "no-Cas strains")) +
  scale_y_continuous(breaks = seq(0, max(data$nplasmids), by = 1), limits = c(0, max(data2$nplasmids))) +
  theme(legend.key.size = unit(1, "cm"), text = element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = xmm, vjust = -5,
               aes(group=set), position = position_dodge(0.6)) +
  errorbar + meanbar
p4

p5 <-
  ggplot(data3, aes(set, nplasmids)) + 
  geom_quasirandom(aes(color=crispr), alpha=1, size=point_size) +
  scale_color_gradient(name = "CRISPR arrays (plasmids)", low = color1_1, high = "red", breaks = seq(0, max(data2$crispr), 1), 
                       limits = c(0,max(data2$crispr))) +
  stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = mm, geom = 'crossbar', color = 'black', 
               linetype = "dashed") +
  theme_classic()  +
  xlab("all strains (Cas9-like)") +
  ylab("Number of plasmids") +
  scale_x_discrete(label=c("Cas9-like strains", "no-Cas9-like strains")) +
  scale_y_continuous(breaks = seq(0, max(data$nplasmids), by = 1), limits = c(0,max(data$nplasmids))) +
  theme(legend.key.height = unit(0.4, "cm"), legend.key.width = unit(1.1, "cm"),
        text = element_text(size=12), legend.position = "none") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = xmm, vjust = -5,
               aes(group=set), position = position_dodge(0.6)) +
  errorbar + meanbar
p5

p1 <-
  ggplot(data, aes(set, nplasmids)) + 
  geom_quasirandom(aes(color=prokka_crispr), alpha=1, size=point_size) +
  scale_color_gradient2(name = "Cas genes (genome)", low = color2_1, high = color2_2, breaks = seq(0, max(data$prokka_crispr), 2), 
                        limits = c(0,max(data$prokka_crispr)), midpoint=5, mid="#537EDB") +
  new_scale_color() +
  geom_quasirandom(aes(color=prokka_arrays), alpha=0.05, size=point_size) +
  scale_color_gradient2(name = "CRISPR arrays (genome)", low = color1_1, high = "darkred", breaks = seq(0, max(data$prokka_arrays), 2),
                       limits = c(0,max(data$prokka_arrays)), midpoint = 5, mid = color1_2) +
  stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = mm, geom = 'crossbar', color = 'black', linetype = "dashed") +
  theme_classic()  +
  xlab("all strains") +
  ylab("Number of plasmids") +
  scale_y_continuous(breaks = seq(0, max(data$nplasmids), by = 1), limits = c(0,max(data$nplasmids))) +
  theme(legend.key.height = unit(0.4, "cm"), legend.key.width = unit(1.1, "cm"),
        text = element_text(size=12), legend.position = "none") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = xmm, vjust = -5,
               aes(group=set), position = position_dodge(0.6)) +
  errorbar + meanbar
p1

p2 <-
  data %>%
  filter (prokka_arrays > 0) %>%
  ggplot(aes(set, nplasmids)) + 
  geom_quasirandom(aes(color=prokka_crispr), alpha=1, size=point_size) +
  scale_color_gradient2(name = "Cas genes (genome)", low = color2_1, high = color2_2, breaks = seq(0, max(data$prokka_crispr), 2), 
                       limits = c(0,max(data$prokka_crispr)), #guide = guide_colourbar(barwidth=1, barheight=8),
                       midpoint=6, mid="#537EDB") +
  new_scale_color() +
  geom_quasirandom(aes(color=prokka_arrays), alpha=0.05, size=point_size) +
  scale_color_gradient2(name = "CRISPR arrays (genome)", low = color1_1, high = "darkred", breaks = seq(0, max(data$prokka_arrays), 2),
                       limits = c(0,max(data$prokka_arrays)), guide = guide_colourbar(barwidth=1, barheight=8),
                       midpoint = 6, mid = color1_2) +
  stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = mm, geom = 'crossbar', color = 'black', linetype = "dashed") +
  theme_classic()  +
  xlab("CRISPR arrays > 0") +
  ylab("Number of plasmids") +
  scale_y_continuous(breaks = seq(0, max(data$nplasmids), by = 1), limits = c(0,max(data$nplasmids))) +
  theme(text = element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank(), legend.position = "none") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = xmm, vjust = -5,
               aes(group=set), position = position_dodge(0.6)) +
  errorbar + meanbar
p2

p2b <-
  data %>%
  filter (prokka_crispr > 0) %>%
  ggplot(aes(set, nplasmids)) + 
  geom_quasirandom(aes(color=prokka_crispr), alpha=1, size=point_size) +
  scale_color_gradient2(name = "Cas genes (genome)", low = color2_1, high = color2_2, breaks = seq(0, max(data$prokka_crispr), 2), 
                        limits = c(0,max(data$prokka_crispr)), #guide = guide_colourbar(barwidth=1, barheight=8),
                        midpoint=6, mid="#537EDB") +
  new_scale_color() +
  geom_quasirandom(aes(color=prokka_arrays), alpha=0.05, size=point_size) +
  scale_color_gradient2(name = "CRISPR arrays (genome)", low = color1_1, high = "darkred", breaks = seq(0, max(data$prokka_arrays), 2),
                        limits = c(0,max(data$prokka_arrays)), guide = guide_colourbar(barwidth=1, barheight=8),
                        midpoint = 6, mid = color1_2) +
  stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = mm, geom = 'crossbar', color = 'black', linetype = "dashed") +
  theme_classic()  +
  xlab("Cas genes > 0") +
  ylab("Number of plasmids") +
  scale_y_continuous(breaks = seq(0, max(data$nplasmids), by = 2), limits = c(0,max(data$nplasmids))) +
  theme(text = element_text(size=12),
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.line.y=element_blank()) +
  stat_summary(fun.data = n_fun, geom = "text", hjust = xmm, vjust = -5,
               aes(group=set), position = position_dodge(0.6)) +
  errorbar + meanbar
p2b

p2c <-
  ggplot(data3, aes(set, nplasmids)) + 
  geom_quasirandom(aes(color=prokka_crispr), alpha=1, size=point_size) +
  geom_quasirandom(aes(color=prokka_arrays), alpha=0.05, size=point_size) +
  scale_color_gradient2(name = "CRISPR arrays (genome)", low = color1_1, high = "darkred", breaks = seq(0, max(data$prokka_arrays), 2),
                        limits = c(0,max(data$prokka_arrays)), midpoint = 5, mid = color1_2) +
  stat_summary(aes(ymin = ..y.., max = ..y..), fun.y = mm, geom = 'crossbar', color = 'black', linetype = "dashed") +
  theme_classic()  +
  xlab("all strains (Cas9-like)") +
  ylab("Number of plasmids") +
  scale_x_discrete(label=c("Cas9 strains-like", "no-Cas9-like strains")) +
  scale_y_continuous(breaks = seq(0, max(data$nplasmids), by = 1), limits = c(0,max(data$nplasmids))) +
  theme(legend.key.height = unit(0.4, "cm"), legend.key.width = unit(1.1, "cm"),
        text = element_text(size=12), legend.position = "none") +
  stat_summary(fun.data = n_fun, geom = "text", hjust = xmm, vjust = -5,
               aes(group=set), position = position_dodge(0.6)) +
  errorbar + meanbar
p2c

pdf("../figures/Fig5.pdf", width=10, height=8, paper='special', onefile=FALSE)
empty <- ggplot() + theme(axis.line = element_blank(), panel.background = element_blank())
figa <- ggarrange(p3, p4, common.legend = TRUE, legend = "top", nrow = 1)
figb <- ggarrange(p1, p2, p2b, common.legend = TRUE, legend = "top", nrow = 1)
figc <- ggarrange(p5, p2c, common.legend = FALSE, legend = "top", nrow = 1)
ggarrange(figa, empty, figb, empty, figc, labels = c("a)", "", "b)", "", "c)"), ncol = 1, heights = c(1, 0.2, 1, 0.2, 1))
dev.off()


