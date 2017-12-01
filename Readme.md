# Multi-transcripts toolbox

## Description

Ce travail a été consacré à la simulation des données RNA-seq et microarrays pour appliquer les normalisations et les standardisations à ces données simulées afin de pouvoir les combiner. 

Il contient 3 modules:

    --microarray_simul : Permet de simuler les  données microarrays  à partir d’un modèle prédéfini. Les données simulées ont un comportement similaire aux données microarrays produites par la plateforme « Affimetrix ».
    -- rnaseq_simul : Permet de simuler les données de comptages RNAseq puis les normaliser. Nous proposons trois méthodes de normalisation des données RNAseq : DESeq2, edgeR , et VOOM.
    -- naseq_microarray_fusion : Permet de standardiser puis fusionner les deux  matrices des données simulées.  Nous proposont également trois méthodes de standardisation : Zscore, Zscore Robuste et la QN. 


## Motivations
Plusieurs approches permettent de mesurer l’expression génique. Il y a la technologie des puces à ADN (microarrays) qui reste jusqu’à aujourd’hui la plus utilisée d’entre elles. Et le séquençage d’ARN qui devient la technologie de choix pour les nouvelles expériences. Cependant, la structure de données et les distributions entre les plates-formes diffèrent. 
Notre objectif est de combiner les deux technologies afin de réaliser des analyses sur les données fusionnées. Pour cela...

...

## Usage in command line
...

## Usage in R scripts
### 1:	microarray_simul.r
Cette fonction permet de simuler des données microarrays à partir d’un modèle prédéfini. Les données simulées ont un comportement similaire aux données microarrays produites par la plateforme « Affimetrix ».
Pour cela, l’utilisateur doit fournir un ensemble des paramètres, ou utiliser ceux disponibles par défaut.
#### Usage: 

" microarray_simul.r" --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 –m1 1.4 --m2 0.8
#### Description des arguments

•	"-gn" ou  "--gene_number": Un nombre entier naturel indiquant  le nombre de gènes dans les données simulées. La valeur par défaut est : --gene_number=10,000
•	"-sn1" ou  "--samples_n1" : Un nombre entier naturel indiquant le nombre d’échantillons  du phénotype 1 (condition 1). La valeur par défaut est : --samples_n1=75
•	
•	"-sn2" ou  "--samples_n2" : Un nombre entier naturel indiquant le nombre d’échantillons  du phénotype 1 (condition 1). La valeur par défaut est : --samples_n2=75
•   "-diff "ou  " --diff_genes_ratio" : Un nombre décimal indiquant le the pourcentage des gènes différentiellement exprimés. Sa valeur par défaut est : --diff_genes_ratio=0.1
•	"-up", ou "--up_ratio" : Un nombre décimal indiquant le pourcentage de gènes surexprimés. Sa valeur par défaut est : --up_ratio=0.5
•	"-m1" ou "--m1": Un nombre décimal  correspondant à la différence moyenne entre la moyenne totale et la moyenne  des gènes différentiellement exprimé avec  des valeurs élevées. Sa valeur par défaut est : --m1=1.4
•	"-m2" ou  "--m2": Un nombre décimal  correspondant à la différence moyenne entre la moyenne totale et la moyenne  des gènes différentiellement exprimé avec  des valeurs peu elevées. Sa valeur par défaut est : --m1=0.8
•	"-s" ou "--seed": Un entier utilisé pour générer un nombre aléatoire par l'ordinateur dans le but de rendre la simulation reproductible
#### Plus de détails 
Si l’utilisateur fournit un nombre décimal au lieu d’un nombre entier  pour les trois premiers paramètres, la valeur sera arrondie. 
La fonction ne sera pas exécutée et retournera un message d’erreur dans les cas suivant :
•	Si un des paramètres numérique ne l’est pas
•	Si un des paramètres numériques est négatif 
•	Si le nombre de gènes à simuler est nul
•	Si le nombre d’échantillons est nul 
•	Si les deux paramètres de proportions ne sont pas compris entre 0 et 1
#### Sortie 
La fonction renvoie une matrice de données avec respectivement  le nombre de lignes et de colonnes spécifié par les paramètres d'entrée "--gene_number et "--samples_n1 + "--samples_n2. 


...

## Deployment and usage in Galaxy workflows
...

## Software dependencies
