

args <- commandArgs(trailingOnly=TRUE)
dossier_data = args[1]
Dossier_Normalisation = args[2]
Dossier_output = args[3]
prefix= args[4]

############################ Cross plateformes

suppressMessages(source("package_loader.R"));
load_it( c("MASS", "DESeq","DESeq2" ,"edgeR","limma","preprocessCore","clusterSim"));
source("normalize_data.R")

 
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


Cross_platformes_normalisation=function(Microarray,RNA,Norm,St,Mn1,Mn2,Rn1,Rn2,common_sample=c("FALSE","TRUE")){

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

	## Cas ou les données sont appariées
	if(common_sample=="TRUE"){
	Cond1_common_samples<-intersect(colnames(RNA[,1:Rn1]),colnames(Microarray[,1:Mn1])) 
	Cond2_common_samples<-intersect(colnames(RNA[,(Rn1+1):ncol(RNA)]),colnames(Microarray[,(Mn1+1):ncol(Microarray)]))
	common_samples=c(Cond1_common_samples,Cond2_common_samples)
	counts=RNA[genes_commun,common_samples]
	Microarray=Microarray[genes_commun,common_samples]
	Mn1<-Rn1<-length(Cond1_common_samples)
	Mn2<-Rn2<-length(Cond2_common_samples)
	} 
	
	
		## Choix de la normalisation des données RNA seq
		R_Norm<-Normalization(counts,Rn1,Rn2,Norm)
		Design=R_Norm[[1]]
		message(" Voulez vous enregistrer la matrice des données RNAseq normalisée ?")
		reponse<-select.list(c("Oui","Non"))
		if( reponse =="Oui"){  save(R_Norm,file=paste0(Dossier_Normalisation,paste0("RNA_",prefix,"_",Norm,"normalisation.RData")))}
		## Récupération de la matrice des comptages normalisée et standarisée
		R<-Standarisation(microarrays,R_Norm$normalized_data,St)
		## Récupération de la matrice des ratio  log intensité standarisée
		M=Standarisation(microarrays,microarrays,St)
		
		# Fusion
		Fusion<-data.frame(M[,1:Mn1],R[,1:Rn1],M[,(Mn1+1):ncol(M)],R[,(Rn1+1):ncol(R)],stringsAsFactors=FALSE) #		   
		colnames(Fusion)=c(paste("Cond1",1:(Mn1+Rn1)),paste("Cond2",1:(Mn2+Rn2)))
		DESIGN=factor(c(rep(1,(Mn1+Rn1)),rep(2,(Mn2+Rn2))))
			   
				return(list(Fusion=Fusion,Design=DESIGN))
           
}       
 
 RNA<-read.table(file=paste0(dossier_data,"RNA_count_",prefix,"_pDE=",1-pi0,".txt"),sep="\t")
Microarray<-read.table(file=paste0(dossier_data,"Array_",prefix,"_pDE=",1-pi0,".txt"),sep="\t")

Cross_platformes<-Cross_platformes_normalisation(Microarray,RNA,Norm,St,n1,n2,n1,n2,common_sample)
write.table( Cross_platformes$Fusion,file=paste0(Dossier_output,paste0("Fusion_",prefix,"_",Norm,"_",St,".txt")),sep="\t",row.names=T, col.names=T) 
