#########################
# Alejandro Rubio, 2019 #
# 16S phylogeny         #
#########################

library(ggplot2)
library(ggtree)

####lectura del árbol con la libreria ggtree

tree<-read.tree("outtree_def") 

####Carga el heatmap, un archivo tsv con las características de cada especie que se quiere incluir en la filogenia

heatmap_tabla<- read.table("heatmap_def_16S.tsv", sep = "\t", header = T, quote = "", 
                           row.names = 1, 
                           stringsAsFactors = FALSE)

####Monta la filogenia inicial

phylo<- ggtree(tree,layout = 'circular', branch.length='none', size= 0.1) + geom_tiplab2(size=0.1) + geom_text2(aes(label=label, subset=!isTip, hjust=1)) #geom_treescale(x=10, y=10) +

####Indica el tamaño del nombre de cada cepa y el ángulo

phylo<- phylo + geom_tiplab(size=1.5, aes(angle=angle))

####Vector con el número de los nodo con un bootstrap inferior a 70

clades<-c(188,187,185,183,173, 158, 157, 156, 154, 139, 136)

####Se incorporan a la filogenia

phylo<- phylo + geom_point2(aes(subset=(node %in% clades)),color="red",size=2)

#####Se monta un vector con las caracteristicas que se van a comparar

q1 <- c(">1000","15-1000","5-15","5","1",
        "97%","94-90%","90-85%", "85-80%", "<80%",
        "A.baumannii","ACB","Other")

#####Se añaden los colores para cada atributo

c1 <- c("#028123", "#38BB5A", "#7AD993", "#B6E7C3", "#DDF4E3",
        "#6F1010","#9B3C3C","#B55E5E","#D69898","#EAC9C9","dodgerblue4","#71B4EA","gray86")

#####Se asocian ambas

names(c1) <- q1

#####Se crea un archivo pdf para guardar la figura, con width y height se puede jugar con el tamaño de la filogenia y leyenda

pdf("../figures/Fig1.pdf", width=10, height=10, paper='special')

#####Se une el heatmap
gheatmap(phylo, heatmap_tabla, offset = 0.5,  width=1, font.size=0,colnames_angle=90, hjust=0.1) + 
  
####Se añade la leyenda
  scale_fill_manual(name= "No.strains     ANI              Species                   ",
                    values=c1, 
                    breaks=q1) + 
  guides(fill=guide_legend(ncol=3)) +
  theme(legend.title = element_text(),
        legend.position="right", 
        legend.key = element_rect(fill="white") ) 

dev.off()

