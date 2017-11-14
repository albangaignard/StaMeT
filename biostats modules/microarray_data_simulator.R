
suppressMessages(source("package_loader.R"))
load_it( c("lattice","MASS","edgeR"))
source("Repertoires.R")
args <- commandArgs(trailingOnly = TRUE)
nGenes=as.numeric(args[1])
n1 = as.numeric(args[2])
n2=as.numeric(args[3])
pi0=as.numeric(args[4])
up=as.numeric(args[5])
muminde1=as.numeric(args[6]);
muminde2=as.numeric(args[7]);


if(is.null(muminde1)){ muminde1=1.4}
if(is.null(muminde2)){ muminde1=0.8}





Simulation.microarray=function(nGenes,n1,n2, ratio=FALSE,pi0, up,  muminde1, muminde2){
## Initialisation des paramètres
		shape2 = 4; lb = 4; ub = 14; lambda1 = 0.13; sdde = 0.5;sdn = 0.4;N=n1+n2+1;lambda2=4; 

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
		## simulation de n valeurs selon la loi bêta, une par gène, d'où découleront ensuite les valeurs pour les différents patients.
        x <- rbeta(nGenes, 2, shape2) 
		## mise à l'échelle (de lb à (lb+ub)): 
		## On transforme ces valeurs pour correspondre aux données réelles variant entre la borne inférieure  lb=4 et la borne supérieure ub=14. x2=lb+ub*x.
		## Ces valeurs représentent le niveau d’expression moyen pour les n gènes.
        x2 <- lb + (ub * x) 
		## initialisation de la matrice avec la bonne taille (ou 1 de trop pour l'éch ref si on est pas en ratio)
    	xdat <- matrix(c(rep(0, nGenes * N)), ncol = N) 
		## "id" des gènes, initialisées à 0 pour tous => pour récupérer info up/down/non de ensuite
		xid <- matrix(c( rep(1, TP_up), rep(-1, TP_down),rep(0, FP)), ncol = 1)		
		LFC <- matrix(c(rep(0, nGenes)), ncol = 1) 
		## on fait une boucle sur les gènes
		for (i in 1:nGenes) {
			alpha <- lambda1 * exp(-lambda1 * x2[i]) 
			## x2[i] contient la valeur "échantillonnée" pour le gène i; on génère aléatoirement N valeurs depuis la loi uniforme (cf algo)
			## Pour chaque valeur de x2, on génère N valeurs uniformément réparties sur un intervalle centré
			xi_val <- runif(N, min = (1 - alpha) * x2[i], max = (1 + alpha) * x2[i]) 
			## Si le gène n'est pas DE
			if (i > TP) { 
            xdat[i, ] <- xi_val 
			LFC[i,]=0
			} else {	## si le gène DE
				## on garde les (n1+1) premières valeurs, qui reflètent les valeurs dans des conditions "normales"/ première condition
				xi1 <- xi_val[1:(n1 + 1)] 
				## on tire au sort la diff de moyenne, selon muminde et lambda2
				 if(( i<(TP_up_up+1)) | ((i>TP_up)&(i<TP_up+TP_down_down+1))) {

				mude <- muminde1 + rexp(1, lambda2) 
				## on tire au sort si gène up ou down, up pour le if, down pour le else
				if (i <= TP_up_up) { 
				fc<-rnorm(n2, mean = mude, sd = sdde)
				LFC[i,]<-mean(fc)
				## on ajoute la partie "up", distribuée selon loi N autour de mude avec toujours la même dispersion sdde
                xi2 <- xi_val[(n1 + 2):N] + fc 
                xid[i] <- 1 # le gène est up
				} else {
					fc<-rnorm(n2, mean = mude, sd = sdde)
					LFC[i,]<-mean(fc)
					## idem que pour les up mais on retranche pour les down
					xi2 <- xi_val[(n1 + 2):N] - fc 
					# gène est down
					xid[i] <- -1 
				}
		
				 xdat[i, ] <- c(xi1, xi2)

			} else {
			
			mude <- muminde2 + rexp(1, lambda2) 
				## on tire au sort si gène up ou down, up pour le if, down pour le else
				if ((i >TP_up_up)&(i<=TP_up)) { 
				fc<-rnorm(n2, mean = mude, sd = sdde)
				LFC[i,]<-mean(fc)
				## on ajoute la partie "up", distribuée selon loi N autour de mude avec toujours la même dispersion sdde
                xi2 <- xi_val[(n1 + 2):N] + fc 
                xid[i] <- 1 # le gène est up
				} else {
					fc<-rnorm(n2, mean = mude, sd = sdde)
					LFC[i,]<-mean(fc)
					## idem que pour les up mais on retranche pour les down
					xi2 <- xi_val[(n1 + 2):N] - fc 
					# gène est down
					xid[i] <- -1 
				}}
	            xdat[i, ] <- c(xi1, xi2)
				}
				}
				
	LFC=LFC*xid

	## annotation des colonnes
    colnames(xdat) <- c("V1", rep("cond1",n1),rep("cond2",n2)) 
	rownames(xdat)=c(paste("Gene.up",1:TP_up),paste("Gene.down",1:TP_down),paste("Gene" ,1:FP))
	## sd de la 1e colonne = "individu" ref pour les ratios
    xsd <- sd(xdat[, 1]) 
	## si sd est strictement positif, on utilise cette valeur pour recréer une matrice de la même taille que xdat, de moyenne 0 et sd=xsd (pour faire le bruit additif)
    if (sdn > 0) { 
        ndat <- matrix(c(rnorm(nGenes * N, mean = 0, sd = sdn)), ncol = N)
		## On ajoute du bruit additif       
	   xdat <- xdat + ndat 
    }
	 ## initialisation de xdata
	xdata <- matrix(c(rep(0, nGenes * (N - 1))), ncol = (N - 1))
	#€ si on est en simulation de ratio, on retranche la ref (correspond à un ratio car on simule données en log)
    if (ratio) { 
        xdata <- xdat[, 2:N] - xdat[, 1]
    } else { # sinon on vire la 1e colonne
        xdata <- xdat[, 2:N]
    }
	
	list(xdata = xdata, xid = xid, xsd = xsd,delta=LFC)
}



Array_data<-Simulation.microarray( nGenes,n1,n2,ratio=FALSE,pi0, up,  muminde1, muminde2)
Microarray<-Array_data$xdata

write.table(Microarray,file=file.path(data_simulation,"Array_simulation.txt"),sep="\t",row.names=T, col.names=T)
