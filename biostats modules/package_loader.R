

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