# Script simulateur RNA-Seq

pack_dispo <- row.names(installed.packages())
if(! "optparse" %in% pack_dispo) install.packages("optparse", repos="https://cloud.r-project.org/")

library(optparse)

option_list = list(
    make_option(c("-gn", "--gene_number"), type="numeric", default=10000, 
              help="Total number of genes that are simulated [default: %default]"),
	make_option(c("-sn1", "--samples_n1"), type="numeric", default=75, 
              help="Number of samples with phenotype 1 [default: %default]"),
	make_option(c("-sn2", "--samples_n2"), type="numeric", default=75, 
              help="Number of samples with phenotype 2 [default: %default]"),
	make_option(c("-up", "--up_ratio"), type="numeric", default=0.5, 
              help="Proportion of up-regulated genes within differentially-expressed genes [default: %default]"),
	make_option(c("-diff", "--diff_genes_ratio"), type="numeric", default=0.1, 
              help="Proportion of differentially expressed genes (genes related to phenotype) within all genes [default: %default]"),
	make_option(c("-fc", "--fc_file"), type="character", default="FC.txt", 
              help="Text file with Fold-Change values for the simulated genes. See documentation for file structure details. [default: %default]"),			  
	make_option(c("-rseq_n", "--rnaseq_norm"), type="character", default="DESeq2", 
              help="Normalisation method for count data, possible values: DESeq2, edgeR, VOOM. [default: %default]"),
    make_option(c("-s", "--seed"), type="numeric", default=NULL, 
              help="Seed value: can be set to make the simulation reproducible. [default: %default]")
)

arg_parser = OptionParser(option_list=option_list)
les_args = parse_args(arg_parser)

if(! "MASS" %in% pack_dispo) install.packages("MASS", repos="https://cloud.r-project.org/")
if(! "edgeR" %in% pack_dispo) install.packages("edgeR", repos="https://cloud.r-project.org/")
library(MASS)
library(edgeR)

counts.simulation <- function(nGenes, n1, n2, pi0, up, fc, seed=NULL){ 
									
	## si un seed a été entré, on le fixe
    if(!is.null(seed)) set.seed(seed)
	
	## moyenne et dispersion
	mu <- runif(nGenes, 500, 3000); disp=rlnorm(nGenes, -1.13, 1)
	replace <- TRUE			
	## Nombre des Faux positifs 																			
	FP <- round(nGenes * pi0) 																			
	## Nombre des vrais positifs																			
	TP <- nGenes - FP 																					
	## types des TP (up & down)																		
	TP_up <- round(TP * up)
	TP_down <- TP - TP_up 

	## vecteur des DE: 1 pour les gènes Up, -1 pour les gènes down et 0 pour les gènes non DE
	DE <- c(rep(1, TP_up), rep(-1, TP_down), rep(0, FP))

    ## LFC, environ moitié positif, moitié négatif (pour les gènes différentiellement exprimés (DE))
	delta <- rep(0, nGenes)
	 	
	# Verifier qu'il n'y a pas de NA 
	# sinon, on sort de la fonction avec message d'erreur
	if(any(is.na(fc))) {message("NA are not allowed in fc_file.") ; return()}

	if((length(fc)==TP) & all(fc[1:TP_up]>1) & all(fc[(TP_up+1):TP]<1)){
     	lfc <- log(fc)
		delta[DE != 0] <- lfc[DE != 0]
	}else{ 
		message("Error: Vector FC is in discordance with parameters pi0 and up:")
		message("---The size of vector FC must be equal to: ", TP)
		message("---The first ", TP_up, " values must be strictly greater than 1")
		message("---The following ", TP_down, " values must be strictly less than 1")
		message("---A default FC file is available when no value is given for fc_file parameter.")
		return()
	}

	 # initialisation
	 h <- rep(TRUE, nGenes) 
	 counts <- matrix(0, nrow = nGenes, ncol = n1+n2) 
	 selected_genes <- true_means <- true_disps <- rep(0, nGenes)
	 left_genes <- 1:length(mu)
	 lambda <- phi <- matrix(0, nrow=nGenes, ncol=n1+n2)
	  
	while(any(h)){
		temp <- sample.int(length(left_genes), sum(h), replace)
		temp <- temp[order(temp)]
		selected_genes[h] <- left_genes[temp]
		if (!replace){
		  left_genes <- left_genes[-temp]
		}
		
		true_means[h] <- mu[selected_genes[h]]
		true_disps[h] <- disp[selected_genes[h]]
		
		lambda[h,] <- matrix(true_means[h], ncol = 1) %*% 
					  matrix(rep(1, n1+n2), nrow = 1) * 
					  cbind(matrix(rep(1, sum(h) * n1), ncol = n1), 
							matrix(rep(exp(delta[h]), n2), ncol = n2))
		
		## moyenne des comptages
		phi[h, ] <- matrix(rep(true_disps[h],n1+n2 ), ncol = n1+n2)
		## dispersion des comptages
		counts[h, ] <- rnegbin(sum(h) * (n1+n2), lambda[h, ], 1/phi[h, ])
	}
	  
	  rownames(counts) <- c(paste("Gene.up", 1:TP_up), paste("Gene.down", 1:TP_down), paste("Gene", 1:FP))
	  delta <- delta/log(2)
	  
	colnames(counts)=c(paste("cond1", 1:n1, sep="_"), paste("cond2", 1:n2, sep="_"))
	list(counts = counts, DE=DE)
}



RNAseq_counts <- counts.simulation(nGenes=les_args$gene_number, n1=les_args$samples_n1, n2=les_args$samples_n2, pi0=1-les_args$diff_genes_ratio, up=les_args$up_ratio, 
                                   fc=read.table(les_args$fc_file, as.is=TRUE, header=TRUE, sep="\t")[, 1, drop=TRUE], seed=les_args$seed)
counts=RNAseq_counts$counts


# NORM
if(! "DESeq2" %in% pack_dispo) install.packages("DESeq2", repos="https://cloud.r-project.org/")
if(! "limma" %in% pack_dispo) install.packages("limma", repos="https://cloud.r-project.org/")
library("DESeq2")
library("limma")
 
Normalization <- function(counts, n1=les_args$samples_n1, n2=les_args$samples_n2, Norm=c("DESeq2", "edgeR", "VOOM")){
    ## Verifier que n1 et n2 soient suprieurs 0
	## si ce n'est pas le cas, on va mettre un design: ~1
    if(n1 & n2){
    	condition <- factor(rep(paste0("Cond", 1:2), c(n1, n2)))
	} else {
	    condition <- rep(1, n1+n2)
	}
	
	# normalisation
	switch(Norm,
	      "DESeq2"={colData <- data.frame(condition, row.names=colnames(counts))
					dds <- DESeqDataSetFromMatrix(counts, colData, design=~condition)
					cds=estimateSizeFactors(dds) 
					#sizeFactors(cds)
					Data_Norm <- counts(cds, normalized=TRUE) 
					Norm_DESeq2_log2=log2(Data_Norm+1)
					output <- Norm_DESeq2_log2},
          "edgeR"={cds <- DGEList(counts, group=condition )
				   cds <- calcNormFactors(cds, method="TMM")
				   lcpm <- cpm(cds, log=TRUE)
				   output <- lcpm},
		  "VOOM"={design=model.matrix(~0+condition)
				  voom_trans<-voom(counts, design,span = 0.5, plot = FALSE,save.plot = FALSE )
				  voom_matrix <-  voom_trans$E
				  output <- voom_matrix})
	return(output)

}

Norm_count <- Normalization(counts, Norm=les_args$rnaseq_norm)

write.table(data.frame(Gene=row.names(Norm_count), Norm_count), file="RNAseq_simulation.txt", sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)


