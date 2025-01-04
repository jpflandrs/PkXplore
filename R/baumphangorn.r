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

# """
# baumphangorn.r is a part of the Web-Site engine PkXplore

# Copyright or © or Copr. UCBL Lyon, France;  
# contributor : [Jean-Pierre Flandrois] ([2024/12/20])
# [JP.flandrois@univ-lyon1.fr]

# This software is a computer program whose purpose is to create the Web-Site PkXplore.

# This software is governed by the [CeCILL|CeCILL-B|CeCILL-C] license under French law and
# abiding by the rules of distribution of free software.  You can  use, 
# modify and/ or redistribute the software under the terms of the [CeCILL|CeCILL-B|CeCILL-C]
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info". 

# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability. 

# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or 
# data to be ensured and,  more generally, to use and operate it in the 
# same conditions as regards security. 

# The fact that you are presently reading this means that you have had
# knowledge of the [CeCILL|CeCILL-B|CeCILL-C] license and that you accept its terms.

# """