#!/usr/bin/env Rscript
library(ape)
library(phangorn)


args <- commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args) < 3) {
  stop("At least 3 arguments must be supplied ",
       call. = FALSE)
}


arbrebrut <- read.tree(args[1])
arbre_avec_racine <- midpoint(arbrebrut, node.labels = "support")
write.tree(arbre_avec_racine, file = args[2], append = FALSE,
           digits = 10, tree.names = FALSE)



farbe <- rep("black", Ntip(arbre_avec_racine))
farbe[grep("QUERY", arbre_avec_racine$tip.label)] <- "red"


#pdf(file = args[2],  width = 27, height = 19, paper = "a4r")
svg(filename=args[3], width = 29, height = 21, pointsize = 32, family = "sans", bg = "white",
    antialias = "default") 

plot(arbre_avec_racine, type = "phylogram",
     no.margin = FALSE,
     show.tip.label = TRUE,
     cex = 0.5,
     plot = TRUE,
     x.lim = NULL,
     y.lim = NULL,
     show.node.label = TRUE,
     tip.color = farbe)

add.scale.bar(cex = 0.8, font = 2, col = "red",lwd = 2, lcol = "red")
dev.off()