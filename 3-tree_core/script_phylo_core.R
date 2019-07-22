#########################
# Alejandro Rubio, 2019 #
# Core phylogeny        #
#########################

library(ggplot2)
library(ggtree)
library(cowplot)

####lectura del árbol con la libreria ggtree

tree<-read.tree("RAxML_bipartitions.Tree_core_def_repl") 

####Carga el heatmap, un archivo tsv con las características de cada especie que se quiere incluir en la filogenia

heatmap_tabla<- read.table("heatmap_def_core.tsv", sep = "\t", header = T, quote = "", 
                           row.names = 1, 
                           stringsAsFactors = FALSE)

####Monta la filogenia inicial

phylo<- ggtree(tree,layout = 'circular', branch.length='none', size= 0.2) + geom_tiplab2(size=0.1)  #geom_treescale(x=10, y=10) +

####Indica el tamaño del nombre de cada cepa y el ángulo

phylo<- phylo + geom_tiplab(size=0.1, aes(angle=angle))


q1 <- c("Blood","Catheter","Inert surface","Osteoarticular","Non-human Host",
        "Other_source","Perianal","Respiratory", "Skin and soft tissue infection", "Urinary and renal fluid", "GROUP_1","GROUP_2",
        "A.baumannii","ACB","OTHER")
#"Non-human host"
#####Se añaden los colores para cada atributo

c1 <- c("red","darkorange","darkturquoise","green","darkred","grey",
        "magenta","darkgreen","blue","yellow3", "orange", "purple","dodgerblue4","lightblue3","gray86")
#"yellowgreen"
#####Se asocian ambas

names(c1) <- q1

pdf("Sup_fig_2.pdf", width=30,height=20, paper='special')

#####Se une el heatmap

phylo2<-gheatmap(phylo, heatmap_tabla, offset = 0.1,  width=0.5, font.size=0 ,colnames_angle=90, hjust=1) +
  
####Se añade la leyenda
  scale_fill_manual(name= "TITULO",
                    values=c1, 
                    breaks=q1) + 
  guides(fill=guide_legend(nrow =10, ncol=3)) +
  theme(legend.text=element_text(size=20),
        legend.title = element_text(size=20, face = "bold"),
        legend.position="right", 
        legend.key = element_rect(fill="white")) 
  #leg1 <- get_legend(phylo2)
  #plot_grid(phylo2, ncol=1, rel_widths=c(1, .1, .1))
phylo2

dev.off()
