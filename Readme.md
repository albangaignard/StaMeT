# Multi-transcripts toolbox

## 1    Description & Motivation

Plusieurs approches permettent de mesurer l’expression génique. Il y a la technologie des puces à ADN (microarrays), qui reste jusqu’à aujourd’hui la plus utilisée d’entre elles, et le séquençage d’ARN qui devient la technologie de choix pour les nouvelles expériences. 

Notre objectif est de combiner les deux technologies afin de réaliser des analyses sur les données fusionnées. Cependant, la nature des données issues de ces deux technologies diffèrent, ce qui qui rend leur combinaison  difficile.

Cet outil permet de simuler des données RNA-seq et microarrays normalisées, de normaliser des données RNA-seq brutes (données de comptage) et de standardiser des données microarrays ou RNA-seq normalisées puis de les fusionner en un tableau unique. 

Il contient 4 modules:

- microarray_simul : Permet de simuler les  données microarrays  à partir d’un modèle prédéfini. Les données simulées ont un comportement similaire aux données microarrays produites par la plateforme « Affymetrix », après normalisation.
- rnaseq_simul : Permet de simuler les données de comptages RNAseq puis de les normaliser. Trois méthodes de normalisation des données RNAseq sont disponibles : DESeq2, edgeR et VOOM.
- normalisation.rna_seq.r: Permet de normaliser les données RNA-seq brutes (données de comptage). Trois méthodes de normalisation des données RNAseq sont disponibles : DESeq2, edgeR et VOOM
- naseq_microarray_fusion : Permet de standardiser puis fusionner des matrices de données RNA-seq ou microarrays normalisées (réelles ou simulées).  Trois méthodes de standardisation sont disponilbles : Zscore, Zscore Robuste et la Quantile Normalisation. 
    


 



## 2	Utilisation sous R 

### 2.1	“microarray_simul.r”

Cette fonction permet de simuler des données microarrays à partir d’un modèle prédéfini. Les données simulées ont un comportement similaire aux données microarrays produites par la plateforme « Affymetrix ».
Pour cela, l’utilisateur doit fournir un ensemble des paramètres, ou utiliser ceux disponibles par défaut.



####    2.1.1	Description des arguments

-	"gene_number": Un nombre entier naturel indiquant le nombre de gènes dans les données simulées. La valeur par défaut est : --gene_number=10,000
-	"samples_n1": Un nombre entier naturel indiquant le nombre d’échantillons du phénotype 1 (condition 1). La valeur par défaut est : --samples_n1=75
-	"samples_n2": Un nombre entier naturel indiquant le nombre d’échantillons du phénotype 2 (condition 2). La valeur par défaut est : --samples_n2=75
-	"diff_genes_ratio": Un nombre décimal indiquant le pourcentage des gènes différentiellement exprimés. Sa valeur par défaut est : --diff_genes_ratio=0.1
-	"up_ratio": Un nombre décimal indiquant le pourcentage de gènes surexprimés. Sa valeur par défaut est : --up_ratio=0.5
-	"m1": Un nombre décimal correspondant à la différence moyenne entre la moyenne totale et la moyenne des gènes différentiellement exprimé avec des valeurs élevées. Sa valeur par défaut est : --m1=1.4
-	"m2": Un nombre décimal correspondant à la différence moyenne entre la moyenne totale et la moyenne des gènes différentiellement exprimé avec des valeurs peu élevées. Sa valeur par défaut est : --m1=0.8
-	"seed": Un entier utilisé pour générer un nombre aléatoire par l'ordinateur dans le but de rendre la simulation reproductible

####    2.1.2	 Plus de détails 

Si l’utilisateur fournit un nombre décimal au lieu d’un nombre entier pour les trois premiers paramètres, la valeur sera arrondie. 

La fonction ne sera pas exécutée et retournera un message d’erreur dans les cas suivant :
-	Si un des paramètres numérique ne l’est pas
-	Si un des paramètres numériques est négatif 
-	Si le nombre de gènes à simuler est nul
-	Si le nombre d’échantillons est nul 
-	Si les deux paramètres de proportions ne sont pas compris entre 0 et 1


####    2.1.3	Sortie 

La fonction renvoie une matrice de données avec respectivement le nombre de lignes et de colonnes spécifié par les paramètres d'entrée "--gene_number et "--samples_n1 + "--samples_n2. 

Les données sont supposées etre semblables aux données microarrays produites par la plateforme « Affymetrix » log2 intensité.

#####   2.1.4	Exécution avec Rscript

RScript " \microarray_simul.r" --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 –m1 1.4 --m2 0.8

###     2.2	“	rnaseq_simul.r”

Cette fonction permet de simuler les données de comptages RNA-seq puis de les normaliser. Trois méthodes de normalisation des données RNA-seq sont disponibles : DESeq2 (1) (2) (3), edgeR (4) (5) et VOOM (6) (7).

####    2.2.1	Description des arguments

-	"gene_number": Un nombre entier naturel indiquant le nombre de gènes dans les données simulées. La valeur par défaut est : --gene_number=10,000
-	"samples_n1" : Un nombre entier naturel indiquant le nombre d’échantillons du phénotype 1 (condition 1). La valeur par défaut est : --samples_n1=75
-	"samples_n2" : Un nombre entier naturel indiquant le nombre d’échantillons du phénotype 2 (condition 2). La valeur par défaut est : --samples_n2=75
-	"diff_genes_ratio" : Un nombre décimal indiquant le pourcentage des gènes différentiellement exprimés. Sa valeur par défaut est : --diff_genes_ratio=0.1
-	"up_ratio" : Un nombre décimal indiquant le pourcentage de gènes surexprimés. Sa valeur par défaut est : --up_ratio=0.5
-	"fc_file" : Un fichier ‘txt’ contenant un vecteur FC (fc= "FC.txt"). Si aucun fichier n’est fourni,Le vecteur des fold-change par défaut sera récupéré, soit tel quel, soit avec un échantiollnage pour avoir les bons nombres des génes DE et up
-	"rnaseq_norm" : Un caractère indiquant la méthode de normalisation des données RNA-seq souhaitée. --rnaseq_norm=’’DESeq2 ‘’ est la valeur par défaut, les alternatives à passer sont ‘’edgeR ‘’ et ‘’VOOM’’. 
-	"seed" : Un entier utilisé pour générer un nombre aléatoire par l'ordinateur dans le but de rendre la simulation reproductible



####    2.2.2	Plus de détails :

A l'instar de la fonction "microarray_simul.r", des vérifications seront faites. Voir paragraphe 2.1.2.

####   2.2.3	Sortie

La fonction renvoie une matrice de données RNA-seq normalisée avec respectivement le nombre de lignes et de colonnes spécifié par les paramètres d'entrée "--gene_number et "--samples_n1 + "--samples_n2. 

####    2.2.4	Exécution avec Rscript

Rscript "\ rnasrq_simul.r" --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 

### 2.3	‘’ normalisation.rna_seq.r’’

Cette fonction permet de normaliser les données RNA-seq brutes (données de comptage). Trois méthodes de normalisation des données RNAseq sont disponibles : DESeq2, edgeR et VOOM

####    2.3.1	Description des arguments

-	"count_file": matrice de comptage des données RNA-seq
-	"design" : un fichier « txt » contenant un vecteur de condition en colonne (variable qualitative) décrivant le plan de l’expérience (condition1/ condition2)
-	"rnaseq_norm" méthode de normalisation des données RNA –seq. --rnaseq_norm=’’DESeq2 ‘’ est la valeur par défaut, les alternatives à passer sont ‘’edgeR ‘’ et ‘’VOOM’’. 

####    2.3.2	Sortie 

La fonction renvoie une matrice de données RNA-seq normalisée suivant la méthode de normalisation choisie

####    2.3.3	Exécution avec Rscript

Rscript -“ \normalisation.rna_seq.r’ " --gene_number --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1

### 2.4	‘’rnaseq_microarray_fusion.r‘’

####    2.4.1	Description des arguments

-	"standardisation" : méthode de standardisation. --standardisation=’’zscore‘’ est la valeur par défaut, les alternatives à passer sont ‘robust_zscore ‘’ et ‘quantil’. 
-	"all_genes" : Un logique (TRUE, FALSE) indiquant si la fonction doit retourner tous les gènes ou bien renvoyer seulement les gènes en commun.--all_genes=TRUE par défaut.
-	"tables" : Un vecteur de caractére contenant les noms des matrices des données à fusionner. Par défaut : "--tables"= " MicroArray_simulation.txt, RNA-seq_simulation.txt"




####    2.4.2	Plus de détails 

Un prétraitement est réalisé afin de rendre les données pertinentes et exploitables comme l’existence des doublons dans les noms de patients entre les tableaux. En effet, S'il y a des doublons, on ajoute un suffixe "_ti" au bout du nom.

Aussi, la vérification de la nature des données. Ils doivent être numériques, dans le cas contraire, ils seront transformés en données numériques. Puis une vérification des données manquantes est effectuée. Dans ce cas-là, la fonction s’arrête et renvoie un message d’erreur.

####    2.4.3	Sortie

La fonction renvoie une matrice de données fusionnées traitées avec la méthode de standardisation choisie.

####      2.4.4	Execution avec Rscript

Rscript “\rnaseq_microarray_fusion.r" "--standardisation ’’zscore‘’ "--tables" " MicroArray_simulation.txt, RNASeq_simulation.txt"

#   3	Bibliographie
1. Huber, Simon Anders and Wolfgang. Differential expression analysis for sequence count data. Genome Biology. 2010.
2. Simon Anders, Wolfgang Huber. Differential expression of RNA-Seq data at the gene level – the DESeq package. Last revision 2016-01-12.
3. DESeq2paper. [http://www-huber.embl.de/DESeq2paper]. [En ligne] 
4. Robinson MD, McCarthy DJ, Smyth GK. edgeR: a Bioconductor package for differential expression analysis of digital gene expression data. Bioinformatics. 2010.
5. edgeR: Empirical Analysis of Digital Gene Expression Data in R. https://bioconductor.org/packages/release/bioc/html/edgeR.html. [En ligne] 
6. Charity W Law, Yunshun Chen, Wei Shi and Gordon K Smyth. voom: precision weights unlock linear model analysis tools for RNA-seq read counts. Genome Biology. 15:R29, 2014.
7. Limma: Linear Models for Microarray Data. http://bioconductor.org/packages/release/bioc/html/limma.html. [En ligne] 




## Deployment and usage in Galaxy workflows
...

## Software dependencies
