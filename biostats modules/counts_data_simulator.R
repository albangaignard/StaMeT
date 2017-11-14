


suppressMessages(source("package_loader.R"))
load_it( c("lattice","MASS","edgeR"))
source("Repertoires.R")
args <- commandArgs(trailingOnly = TRUE)
nGenes=as.numeric(args[1])
n1 = as.numeric(args[2])
n2=as.numeric(args[3])
pi0=as.numeric(args[4])
up=as.numeric(args[5])
vect_fc=args[6]

    seed=50
	set.seed(seed)
	mu=runif(nGenes,500,3000)
	disp=rlnorm(nGenes,-1.13,1)
	replace=TRUE

if ( !is.null(vect_fc)){
fc=read.table(vect_fc)[,1]}


counts.simulation<- function(nGenes, n1,n2,pi0, mu, disp,up,replace ,fc){ 
									
	  ## Nombre des Faux positifs 																			
	  FP <- round(nGenes * pi0) 																			
	  ## Nombre des vrais positifs																			
	  TP <- nGenes - FP 																					
	  ## types des TP ( up & down)																		
	  TP_up <- round(TP * up)
	  TP_down <- TP - TP_up 
	  ## TP avec des LFC élevés en valeur absolue
	  TP_up_up=round((TP_up/2))
	  TP_down_down=round((TP_down/2))
	  #//.......................................................................................................//##
	  ## vecteur des DE: 1 pour les gènes Up, -1 pour les gènes down et 0 pour les gènes non DE
	  de <- c( rep(1, TP_up), rep(-1, TP_down),rep(0, FP))
	   #//.......................................................................................................//##
	    ## LFC, environ moitié positif, moitié négatif (pour les gènes différentiellement exprimés (DE))
	  delta <- rep(0, nGenes)
	  ## choix du FC
	  seed=1234
	  set.seed(seed)
	 if (TP_up_up==0){
		FC=c(sort(rnorm(TP_up,1.4,0.2),decreasing = TRUE),
		sort(rnorm(TP_down,-1.4,0.2),decreasing = FALSE),rep(0,FP))} 
		if(TP_up_up>0){
		FC=c(sort(rnorm(TP_up,1.4,0.2),decreasing = TRUE)[1:TP_up_up],
		sort(rnorm(TP_up,0.6,0.3),decreasing = TRUE)[1:(TP_up-TP_up_up)],
		sort(rnorm(TP_down,-1.4,0.2),decreasing = FALSE)[1:TP_down_down],
		sort(rnorm(TP_down,-0.6,0.3),decreasing = FALSE)[1:(TP_down-TP_down_down)],rep(0,FP))}
		FC=exp(FC*log(2))
	
	  ## Vérification du vecteur FC entré par l'utilisateur

	
	  if(!is.null(fc)){
		
	     if((length(fc)==nGenes)& all(fc[1: TP_up]>1.4) & all(fc[1: TP_up_up]>2.5)& all(fc[(TP_up_up+1):TP_up]<2.5) &  all(fc[(TP_up+1):TP]>0.2) & all(fc[(TP_up+1):TP]<0.85)  & all(fc[(TP_up+1):(TP_up+TP_down_down)]<0.5)){
     		lfc <- log(fc)
		delta[de != 0] <- lfc[de != 0]
	     }else{
	 	message(" Le vecteur FC que vous avez entré ne respecte pas l ordre et la structure des données ou il n'a pas les bonnes dimensions ") 
	  	message(" Nous le remplaçons par un vecteur FC qui remplie ces conditions")
		lfc <- log(FC)
  		delta[de != 0] <- lfc[de!= 0]
		}}else{
 			lfc <- log(FC)
  			delta[de != 0] <- lfc[de!= 0]
	  }

	   
#//.......................................................................................................//##
	  ## h = vector indicating which pseudo-genes to re-simulate,0
	  h <- rep(TRUE, nGenes) 
	  counts <- matrix(0, nrow = nGenes, ncol = n1+n2) 
		  selected_genes <- true_means <- true_disps <- rep(0, nGenes)
	  left_genes <- 1:length(mu)
	  lambda <- phi <- matrix(0, nrow = nGenes, ncol =n1+n2 )
	  
	  while(any(h)){
		temp <- sample.int(length(left_genes), sum(h), replace)
		temp <- temp[order(temp)]
		selected_genes[h] <- left_genes[temp]
		if (replace == FALSE){
		  left_genes <- left_genes[-temp]
		}
		
		true_means[h] <- mu[selected_genes[h]]
		true_disps[h] <- disp[selected_genes[h]]
		
		lambda[h,] <- matrix(true_means[h], ncol = 1) %*% 
					  matrix(rep(1, n1+n2), nrow = 1) * 
					  cbind(matrix(rep(1, sum(h) * n1), ncol = n1), 
							matrix(rep(exp(delta[h]), n2), ncol = n2))
		
		## moyenne des comptages
		
		phi[h,] <- matrix(rep(true_disps[h],n1+n2 ), ncol = n1+n2)
		## disperssion des comptages
		
		counts[h,] <- rnegbin(sum(h) * (n1+n2), lambda[h,], 1 / phi[h,])
		h <- (rowSums(cpm(counts) > 2) < n1)

	  }
	  if(any(rowSums(cpm(counts) > 2) < n1 ))
		print("Erreur: Impossible de simuler des données: certains gènes ne sont pas exprimés.")
	  rownames(counts)=c(paste("Gene.up",1:TP_up),paste("Gene.down",1:TP_down),paste("Gene" ,1:FP))
	  delta <- delta / log(2)
	  
	colnames(counts)=c(rep("cond1",n1),rep("cond2",n2))
	  list(counts = counts,   de = de)
		  
	
}


		 seed=50
		set.seed(seed)
		replace = TRUE
	
		mu=runif(nGenes,500,3000)
		disp=rlnorm(nGenes,-1.13,1)


#



RNAseq_counts<-counts.simulation(nGenes, n1,n2,pi0,up,fc)
counts=RNAseq_counts$counts
write.table(counts, file=file.path(data_simulation,"RNA_count_simulation.txt"),sep="\t",row.names=T, col.names=T)



