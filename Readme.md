# Multi-transcripts toolbox

## Description

Cet outil permet de simuler des données RNA-seq et microarrays normalisées, de normaliser des données RNA-seq brutes (données de comptage) et de standardiser des données microarrays ou RNA-seq normalisées puis de les fusionner en un tableau unique. 

Il contient 4 modules:

    -- microarray_simul : Permet de simuler les  données microarrays  à partir d’un modèle prédéfini. Les données simulées ont un comportement similaire aux données microarrays produites par la plateforme « Affymetrix ».
    -- rnaseq_simul : Permet de simuler les données de comptages RNAseq puis de les normaliser. Trois méthodes de normalisation des données RNAseq sont disponibles : DESeq2, edgeR et VOOM.
    -- normalisation.rna_seq.r: Permet de normaliser les données RNA-seq brutes (données de comptage). Trois méthodes de normalisation des données RNAseq sont disponibles : DESeq2, edgeR et VOOM
    -- naseq_microarray_fusion : Permet de standardiser puis fusionner des matrices de données RNA-seq ou microarrays normalisées (réelles ou simulées).  Trois méthodes de standardisation sont disponilbles : Zscore, Zscore Robuste et la Quantile Normalisation. 
    
    

    


## Motivations

Plusieurs approches permettent de mesurer l’expression génique. Il y a la technologie des puces à ADN (microarrays), qui reste jusqu’à aujourd’hui la plus utilisée d’entre elles, et le séquençage d’ARN qui devient la technologie de choix pour les nouvelles expériences. Cependant, la structure de données et les distributions entre les plates-formes diffèrent. 
Notre objectif est de combiner les deux technologies afin de réaliser des analyses sur les données fusionnées.

 



...

## Usage in command line
...

## Usage in R scripts
### 1:	microarray_simul.r

Cette fonction permet de simuler des données microarrays à partir d’un modèle prédéfini. Les données simulées ont un comportement similaire aux données microarrays produites par la plateforme « Affymetrix ».
Pour cela, l’utilisateur doit fournir un ensemble des paramètres, ou utiliser ceux disponibles par défaut.


#### 1:1: Description des arguments

-	"-gn" ou  "--gene_number":              Un nombre entier naturel indiquant  le nombre de gènes dans les données simulées. La valeur par défaut est : --gene_number=10,000

-	"-sn1" ou  "--samples_n1":              Un nombre entier naturel indiquant le nombre d’échantillons  du phénotype 1 (condition 1). La valeur par défaut est : --samples_n1=75

-	"-sn2" ou  "--samples_n2":              Un nombre entier naturel indiquant le nombre d’échantillons  du phénotype 2 (condition 2). La valeur par défaut est : --samples_n2=75

-   "-diff "ou  "--diff_genes_ratio":       Un nombre décimal indiquant le pourcentage des gènes différentiellement exprimés. Sa valeur par défaut est : --diff_genes_ratio=0.1

-	"-up", ou "--up_ratio":                 Un nombre décimal indiquant le pourcentage de gènes surexprimés. Sa valeur par défaut est : --up_ratio=0.5

-	"-m1" ou "--m1":                        Un nombre décimal  correspondant à la différence moyenne entre la moyenne totale et la moyenne  des gènes différentiellement exprimé avec  des valeurs élevées. Sa valeur par défaut est : --m1=1.4

-	"-m2" ou  "--m2":                       Un nombre décimal  correspondant à la différence moyenne entre la moyenne totale et la moyenne  des gènes différentiellement exprimé avec  des valeurs peu elevées. Sa valeur par défaut est : --m1=0.8

-	"-s" ou "--seed":                       Un entier utilisé pour générer un nombre aléatoire dans le but de rendre la simulation reproductible

#### 1.2: Plus de détails 

- Si l’utilisateur fournit un nombre décimal au lieu d’un nombre entier pour les trois premiers paramètres, la valeur sera arrondie. 
- La fonction ne sera pas exécutée et retournera un message d’erreur dans les cas suivant :

-   ** 	Si un des paramètres numériques ne l’est pas

-   **  Si un des paramètres numériques est négatif   

-   **  Si le nombre de gènes à simuler est nul

-   **  Si le nombre d’échantillons est nul 

-   **  Si les deux paramètres de proportions ne sont pas compris entre 0 et 1
    
#### 1.3: Sortie 

La fonction renvoie une matrice de données avec respectivement  le nombre de lignes et de colonnes spécifié par les paramètres d'entrée "--gene_number et "--samples_n1 + "--samples_n2. 
Les données sont supposées  etre semblables aux données microarrays produites par la plateforme « Affimetrix » log2 intensité.
#### 1.4	Exécution de la fonction avec Rscript

RScript " \microarray_simul.r" --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 –m1 1.4 --m2 0.8

### 2 	rnaseq_simul.r

#### 2.1	Description des arguments
•	"-gn" ou  "--gene_number": Un nombre entier naturel indiquant  le nombre de gènes dans les données simulées. La valeur par défaut est : --gene_number=10,000

•	"-sn1" ou  "--samples_n1" : Un nombre entier naturel indiquant le nombre d’échantillons  du phénotype 1 (condition 1). La valeur par défaut est : --samples_n1=75

•	"-sn2" ou  "--samples_n2" : Un nombre entier naturel indiquant le nombre d’échantillons  du phénotype 1 (condition 1). La valeur par défaut est : --samples_n2=75

•	"-diff "ou  " --diff_genes_ratio" : Un nombre décimal indiquant le the pourcentage des gènes différentiellement exprimés. Sa valeur par défaut est : --diff_genes_ratio=0.1

•	"-up", ou "--up_ratio" : Un nombre décimal indiquant le pourcentage de gènes surexprimés. Sa valeur par défaut est : --up_ratio=0.5

•	"-fc" ou "--fc_file" :   Un fichier ‘txt’ contenant un vecteur FC (fc= "FC.txt"). Si aucun fichier n’est fourni,Le vecteur des fold-change par défaut sera récupéré, soit tel quel, soit avec un échantiollnage pour avoir les bons nombres des génes DE et up

•	"-rseq_n", "--rnaseq_norm" : Un caractère indiquant la méthode de normalisation des données RNA-seq souhaitée.  --rnaseq_norm=’’DESeq2 ‘’ est la valeur par défaut,  les alternatives à passer sont ‘’edgeR ‘’ et ‘’VOOM’’. 

•	"-s", "--seed" : Un entier utilisé pour générer un nombre aléatoire par l'ordinateur dans le but de rendre la simulation reproductible




#### 2.2	Plus de détails :
A l'instar de la fonction "microarray_simul.r", des vérifications seront faites. (Voir le paragraphe 1.2: Plus de détails)

#### 2.3	Sortie
La fonction renvoie une matrice de données RNA-seq normalisées avec respectivement le nombre de lignes et de colonnes spécifié par les paramètres d'entrée "--gene_number et "--samples_n1 + "--samples_n2. 

#### 2.4	Exécution de la fonction avec Rscript

Rscript "\ rnasrq_simul.r" --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 


### 3. normalisation.rna_seq.r

#### 3.1 Description des arguments


•	"-count" ou  "--count_file": matrice de comptage des données RNA-seq

•	"-des" ou "--design" : un fichier « txt » contenant un vecteur de condition en colonne  (variable qualitative) décrivant le plan de l’expérience (condition1/ condition2)

•	"-rseq_n" ou  "--rnaseq_norm"  méthode de normalisation des données RNA –seq. --rnaseq_norm=’’DESeq2 ‘’ est la valeur par défaut,  les alternatives à passer sont ‘’edgeR ‘’ et ‘’VOOM’’. 

#### 3.2	Sortie 

La fonction renvoie une matrice de données RNA-seq normalisée suivant la méthode de normalisation choisie

#### 3.3	Exécution de la fonction avec Rscript

Rscript -“ \normalisation.rna_seq.r’ " --gene_number  --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1

### 4	rnaseq_microarray_fusion.r

#### 4.1	Description des arguments

•	"-stand" ou "--standardisation" : méthode de standardisation. --standardisation=’’zscore‘’ est la valeur par défaut,  les alternatives à passer sont ‘robust_zscore ‘’ et ‘quantil’

•	-all"ou "--all_genes" : Un logique (TRUE, FALSE) indiquant si la fonction doit retourner tous les gènes ou bien renvoyer seulement les gènes en commun.--all_genes=TRUE par défaut.

•	"-t" ou  "--tables" : Un vecteur de caractére contenant les noms des matrices des données à fusionner. Par défaut : "--tables"= " MicroArray_simulation.txt, RNASeq_simulation.txt"


#### 4.2	Plus de détails 

Un prétraitement est réalisé afin de rendre les données pertinentes et exploitables comme l’existence des doublons dans les noms de patients entre les tableaux. En effet, S'il y a des doublons, on ajoute un suffixe "_ti" au bout du nom.

Aussi, la vérification de la nature des données. Ils  doivent être numériques, dans le cas contraire, ils seront transformés en données numériques. Puis une vérification des données manquantes est effectuée. Dans ce cas-là, la fonction s’arrête et renvoie un message d’erreur.

#### 4.3	Sortie

La fonction renvoie une matrice de données  fusionnées traitées avec la méthode de  standardisation choisie.

####4.4	Execution avec Rscript

Rscript  “\rnaseq_microarray_fusion.r " "--standardisation ’’zscore‘’   "--tables"  " MicroArray_simulation.txt, RNASeq_simulation.txt"




## Deployment and usage in Galaxy workflows
...

## Software dependencies
