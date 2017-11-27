# permet à l'utilisateur de charger tous les paquets souhaités, en les installant depuis  CRAN ou Bioconductor 
# Prend comme entrée un vecteur de chaînes de caractères contenant des noms de paquets.


load_it = function(pack, update=FALSE) {
	invisible(sapply(pack, function(package) {
		if(!require(package, quietly=T, character.only=T)){
				if(length(find("biocLite")) < 1) {
						source("http://bioconductor.org/biocLite.R")
				}
				if(update) {
						biocLite(package)	
				} else {
						biocLite(package, suppressUpdates=TRUE)
				}

				library(package, quietly=T, character.only=T)
			}
		}))
}