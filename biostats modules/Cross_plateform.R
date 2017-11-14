############################ Cross plateformes

suppressMessages(source("package_loader.R"))
load_it( c("MASS", "DESeq","DESeq2" ,"edgeR","limma","preprocessCore","clusterSim"))
source("Directories.R")
source("normalize_data.R")


args <- commandArgs(trailingOnly=TRUE)

prefix=args[1]
Norm=args[2]
St=args[3]
Mn1=as.numeric(args[4])
Mn2=as.numeric(args[5])
Rn1=as.numeric(args[6])
Rn2=as.numeric(args[7])
 
Standarisation=function(Ref,dat,St=c("Zscor","med_centering","Robust_Zscor","QN")){
dat<-as.matrix(dat)
if(St=="Zscor"){
	##################### Z scor #########################
	#((x-mean)/sd)
	Zscore=data.Normalization (dat,type="n1",normalization="row")
	output=Zscore} 
		if( St=="med_centering"){
		##################### median centering #############
		centrage_Med<-sapply(1:nrow(dat),function(i){(dat[i,] - median(dat[i,]))})
		centrage_Med=t(	centrage_Med)
		rownames(centrage_Med)=rownames(dat)
		output=centrage_Med} 
			if(St=="Robust_Zscor"){
			##################### Robust Z scor #####################
			#type=n2 : (x-median)/mad)
			Zscore_Robust<-data.Normalization (dat,type="n2",normalization="row")
			output=Zscore_Robust} else {
				#  quantile Normalisation
				# Create a quantile normalized dataset.
				# Create a target object for the quantile normalization.
				ref= data.frame(t(Ref))
				target = data.frame(t(dat))
				targ = normalize.quantiles.determine.target(    data.matrix(ref),    target.length=nrow(ref))
				# Quantile normalize the data, against the reference distribution,
				qn = data.matrix(normalize.quantiles.use.target( data.matrix(target),targ,copy=F))
				QN<-t(qn)
				rownames(	QN)=rownames(dat)
				output=QN}
	return(output)
}


Cross_platformes_normalisation=function(Microarray,RNA,Norm,St,Mn1,Mn2,Rn1,Rn2){

	##  vérifier et supprimer les valeurs manquantes NA
	 Microarray<-Microarray[apply(Microarray, 1, function(x) !any(is.na(x))),]
	RNA<-RNA[apply(RNA, 1, function(x) !any(is.na(x))),]
	## Enlever les gènes qui ont une proportion de zéros supérieure à 30%
	proportion<-apply(RNA, 1, function(x)( sum(x==0)/length(x))*100)
	RNA<-RNA[proportion<30,]

	proportion<-apply(Microarray, 1, function(x)( sum(x==0)/length(x))*100)
	Microarray<-Microarray[proportion<30,]
	## Détecter les gènes communs: filtrer les gènes dans les données de comptage pour inclure uniquement ceux dans les données de microarray
	genes_RNA=rownames(RNA)
	genes_microarray=rownames(Microarray)
	genes_commun<-intersect(genes_RNA,genes_microarray)

	counts=RNA[genes_commun,]
	microarrays=Microarray[genes_commun,]

	
	
	
		## Choix de la normalisation des données RNA seq
		R_Norm<-Normalization(counts,Rn1,Rn2,Norm)
		Design=R_Norm[[1]]
		save(R_Norm,file=paste0(data_normalisation,paste0("RNA_",prefix,"_",Norm,"normalisation.RData")))
		## Récupération de la matrice des comptages normalisée et standarisée
		R<-Standarisation(microarrays,R_Norm$normalized_data,St)
		## Récupération de la matrice des ratio  log intensité standarisée
		M=Standarisation(microarrays,microarrays,St)
		
		# Fusion
		Fusion<-cbind(M[,1:Mn1],R[,1:Rn1],M[,(Mn1+1):ncol(M)],R[,(Rn1+1):ncol(R)]) #		   
		colnames(Fusion)=c(paste("Cond1",1:(Mn1+Rn1)),paste("Cond2",1:(Mn2+Rn2)))
		DESIGN=factor(c(rep(1,(Mn1+Rn1)),rep(2,(Mn2+Rn2))))
			   
				return(list(Fusion=Fusion,Design=DESIGN))
           
}       

 
 RNA<-read.table(file=file.path(data_simulation,paste0("RNA_count_",prefix,".txt")),sep="\t")
Microarray<-read.table(file=file.path(data_simulation,paste0("Array_",prefix,".txt")),sep="\t")

Cross_platformes<-Cross_platformes_normalisation(Microarray,RNA,Norm,St,Mn1,Mn2,Rn1,Rn2)

write.table( Cross_platformes$Fusion,file=file.path(output,paste0("Fusion_",prefix,"_",Norm,"_",St,".txt")),sep="\t",row.names=T, col.names=T) 
