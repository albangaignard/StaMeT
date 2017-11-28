# Script fusion données microarray / RNA-Seq

pack_dispo <- row.names(installed.packages())
if(! "optparse" %in% pack_dispo) install.packages("optparse", repos="https://cloud.r-project.org/")

library(optparse)

option_list = list(
    make_option(c("-stand", "--standardisation"), type="character", default="zscore", 
              help="Standardisation method, possible values: zscore, quantile, robust_zscore [default: %default]"),
	make_option(c("-all", "--all_genes"), type="logical", default=TRUE, 
              help="Should all genes be returned? FALSE will return only common genes between tables [default: %default]"),
	make_option(c("-t", "--tables"), type="character", default="MicroArray_simulation.txt;RNASeq_simulation.txt", 
              help="Path to the text files with data to merge, separated by semi-colons. See documentation for details. [default: %default]")		
)

arg_parser = OptionParser(option_list=option_list)
les_args = parse_args(arg_parser)




if(!require(“edgeR”, quietly=T, character.only=T)){
    source("http://bioconductor.org/biocLite.R")                  
    biocLite("edgeR")           
    library(“edgeR”, quietly=T, character.only=T)
}

if(!require(“DESeq2”, quietly=T, character.only=T)){
    source("http://bioconductor.org/biocLite.R")                  
    biocLite("DESeq2")           
    library(“DESeq2”, quietly=T, character.only=T)

}


if(! "preprocessCore" %in% pack_dispo) install.packages("preprocessCore", repos="https://cloud.r-project.org/")
if(! "clusterSim" %in% pack_dispo) install.packages("clusterSim", repos="https://cloud.r-project.org/")
library("MASS")


library("preprocessCore")
library("clusterSim")


meth <- les_args$standardisation
nom_tables <- strsplit(les_args$tables, ";")[[1]]

### VERIFICATION DES PARAMETRES ENTRES PAR L'UTILISATEUR


###

# import des tableaux
l_tables <- lapply(nom_tables, function(nom_fich) read.table(nom_fich, sep="\t", header=TRUE, as.is=TRUE, row.names=1))
names(l_tables) <- sub("\\..+$", "", basename(nom_tables))
nb_t <- length(l_tables)

# Fonction standardisation 
Standardisation=function(dat, St=c("zscore", "robust_zscore", "quantile"), Ref=l_tables[[1]]){
	
	dat <- as.matrix(dat)
	
	switch(St, 
           "zscore"={output <- data.Normalization(dat, type="n1", normalization="row")},
		   "robust_zscore"={output <- data.Normalization(dat, type="n2", normalization="row")},
		   "quantile"={ref <- data.frame(t(Ref))
				       target <- data.frame(t(dat))
				       targ <- normalize.quantiles.determine.target(data.matrix(ref), target.length=nrow(ref))
				       # Quantile normalize the data, against the reference distribution
				       qn = data.matrix(normalize.quantiles.use.target(data.matrix(target), targ, copy=FALSE))
				       output <- t(qn)
				       rownames(output) <- rownames(dat)})
	return(output)
}

# S'il y a des doublons dans les noms de patients entre les tableaux, on ajoute "_ti" au bout du nom de colonne (i étant le numéro du tableau dans la liste)
# on récupère aussi le nom des éventuels tableaux avec des patients en double
cnames <- lapply(l_tables, colnames)
t_noms_tables <- table(stack(cnames))
t_pb_dupes <- names(cnames)[colSums(t_noms_tables>1)>0] # !!!! VOIR COMMMENT UTILISER L'INFORMATION !!!!
if(any(rowSums(t_noms_tables>0)>1)) {
	l_tables <- setNames(lapply(seq_along(l_tables), 
	                     function(numt) {
						    tab <- l_tables[[numt]]
							colnames(tab) <- paste0(colnames(tab), "_t", numt)
							tab
						 }), 
	                     names(l_tables))
}	

# Standardisation puis fusion de tous les tableaux
Fusion <- Reduce(function(tab1, tab2) {
					fus <- merge(tab1, tab2, by="row.names", all=les_args$all_genes)
					row.names(fus) <- fus[, 1]
					fus[, 1] <- NULL
					return(fus)
				 }, lapply(l_tables, Standardisation, St=meth))

write.table(data.frame(Gene=row.names(Fusion), Fusion), file="Fusion.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE) 

# Remplacement des valeurs manquantes par la médiane pour l'ACP (summary plot), avec transposition du tableau au préalable
tFus <- sapply(as.data.frame(t(Fusion)), function(x) {x[is.na(x)] <- median(x, na.rm=TRUE); x})

# ACP puis plot PP1
acp <- prcomp(tFus, center=TRUE, scale=TRUE)
png("Summary_Plot_Fusion.png")
 par(mar=c(10, 4, 4, 2), xpd=TRUE)
 plot(acp$x[, 1:2], pch=19, las=1, col=rep(rainbow(nb_t), sapply(l_tables, ncol)), main = "Fusion Summary Plot (PCA)")
 par(plt=c(0.15, 0.3, 0.15, 0.18), new=TRUE)
 plot(0:1, 0:1, type="n", axes=FALSE, xlab="", ylab="")
 legend("topleft", bty="n", col=rainbow(nb_t), legend=names(l_tables), pch=19, ncol=nb_t%/%3+sign(nb_t%%3))
dev.off()







