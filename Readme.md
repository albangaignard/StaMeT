# Multi-transcripts toolbox

## I.	Description and Motivation

Several approaches exist for studying and measuring gene expression. Microarrays, which is still the most used, and RNA sequencing, which are gaining ground and becoming the technology of choice for new experiments.
The aim is to perform statistical analysis on merged data from these two technologies, despite their differences in structure and nature. 

This tool can simulate RNA-seq and microarrays normalized data, normalize raw RNA-seq data (count data), and standardize microarrays and/or RNA-seq data and then merge into a single data set.
It contains 4 modules:

-	microarray_simul.r: it allows to generate microarrays data. The simulated data have similar characteristics compared to the microarrays data produced by the Affymetrix® platform, after normalization;
-	rnaseq_simul.r: it allows to simulate the RNA-seq count data and then normalizes them. There are three normalization methods implemented: DESeq2, edgeR, and VOOM;
-	normalisation.rna_seq.r: it allows to normalize RNA-seq count data. Three normalization methods are available: DESeq2, edgeR and VOOM;
-	rnaseq_microarray_fusion.r: it allows to standardize and then merge standardized RNA-seq and/or microarrays data (real or simulated). Three standardization methods are available: Z-score, Z-score Robust and Quantile Normalization.


## II.	Deployment and usage in R 

### II.1	 “microarray_simul.r”

This function allows to generate microarrays data with two conditions and known characteristics. These data have similar behavior as those obtained with Affymetrix® platform, after normalization (Log2 intensity).



####    II.1.1	 Arguments
-	gene_number: an integer specifying the number of genes to be simulated. Default to 10,000;
-	samples_n1: an integer specifying the number of phenotype 1 (condition 1) samples to be simulated. Default to 75;
-	samples_n2: an integer specifying the number of phenotype 2 (condition 2) samples to be simulated. Default to 75;
-	diff_genes_ratio: the proportion of differentially expressed genes. Default to 0.1;
-	up_ratio: the proportion of up-regulated genes among differentially expressed genes. Default to 0.5;
-	m1: a decimal number corresponding to average expression difference between condition 2 and condition 1. Default to 1.4;
-	m2: similar to m1, it allows to have 2 levels of difference, for example high and moderate. Default to 0.8;
-	seed: an integer used as seed for generating random number, it permits to generate reproducible data. By default, none is set.

####    II.1.2	Note and details

If the user supplies a decimal number instead of an integer for the first three parameters, the value will be rounded.
The function will not run and will return error messages in the following cases:

-	one of the numeric parameters is not numeric;
-	one of the integer parameters is negative; 
-	the number of genes to be simulated is zero;
-	the number of samples is zero;
-	at least one proportion parameter is not between 0 and 1.



####    II.1.3	Value
The output is a tab-delimited text file containing a dataset with gene_number rows and (samples_n1+samples_n2+1) columns. The first column contains gene names. First ones are the up–regulated genes, then down-regulated genes, then the genes that are not differentially expressed. 


#####   II.1.4	Usage in Rscript

RScript microarray_simul.r --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 –m1 1.4 --m2 0.8 --seed 123

```

###     2.2	“	rnaseq_simul.r”

Cette fonction permet de simuler les données de comptages RNA-seq puis de les normaliser. Trois méthodes de normalisation des données RNA-seq sont disponibles : DESeq2 (1) (2) (3), edgeR (4) (5) et VOOM (6) (7).

####    2.2.1	Description des arguments

-	`gene_number`: Un nombre entier naturel indiquant le nombre de gènes dans les données simulées. La valeur par défaut est 10000
-	`samples_n1` : Un nombre entier naturel indiquant le nombre d’échantillons du phénotype 1 (condition 1). La valeur par défaut est 75
-	`samples_n2` : Un nombre entier naturel indiquant le nombre d’échantillons du phénotype 2 (condition 2). La valeur par défaut est 75
-	`diff_genes_ratio` : Un nombre décimal indiquant le pourcentage des gènes différentiellement exprimés. Sa valeur par défaut est 0.1
-	`up_ratio` : Un nombre décimal indiquant le pourcentage de gènes surexprimés. Sa valeur par défaut est 0.5
-	`fc_file` : Le chemin vers un fichier contenant un vecteur de fold-changes en colonne, avec une en-tête. Si aucun fichier n’est fourni, le vecteur des fold-changes par défaut sera pris, soit tel quel, soit avec un échantillonnage pour avoir les bons nombres des gènes DE et up
-	`rnaseq_norm` : Un caractère indiquant la méthode de normalisation des données RNA-seq souhaitée. "DESeq2" est la valeur par défaut, les alternatives "edgeR" et "VOOM". 
-	`seed` : Un entier utilisé pour générer un nombre aléatoire par l'ordinateur dans le but de rendre la simulation reproductible



####    2.2.2	Plus de détails :

A l'instar de la fonction "microarray_simul.r", des vérifications seront faites. Voir paragraphe 2.1.2.

####   2.2.3	Sortie

La fonction renvoie une matrice de données RNA-seq normalisée avec respectivement le nombre de lignes et de colonnes spécifié par les paramètres d'entrée `--gene_number` et `--samples_n1` + `--samples_n2`. 

####    2.2.4	Exécution avec Rscript
```
Rscript rnasrq_simul.r --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 
```

### 2.3	‘’ normalisation.rna_seq.r’’

Cette fonction permet de normaliser les données RNA-seq brutes (données de comptage). Trois méthodes de normalisation des données RNAseq sont disponibles : DESeq2, edgeR et VOOM

####    2.3.1	Description des arguments

-	`count_file` : matrice de comptage des données RNA-seq
-	`design` : un fichier « txt » contenant un vecteur de condition en colonne (variable qualitative) décrivant le plan de l’expérience (condition1/ condition2)
-	`rnaseq_norm` : méthode de normalisation des données RNA –seq. "DESeq2" est la valeur par défaut, les alternatives sont "edgeR" et "VOOM". 

####    2.3.2	Sortie 

La fonction renvoie une matrice de données RNA-seq normalisée suivant la méthode de normalisation choisie

####    2.3.3	Exécution avec Rscript
```
Rscript normalisation.rna_seq.r --gene_number --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1
```
### 2.4	‘’rnaseq_microarray_fusion.r‘’

####    2.4.1	Description des arguments

-	`standardisation` : méthode de standardisation. --standardisation=’’zscore‘’ est la valeur par défaut, les alternatives à passer sont ‘robust_zscore ‘’ et ‘quantil’. 
-	`all_genes` : Un logique (TRUE, FALSE) indiquant si la fonction doit retourner tous les gènes ou bien renvoyer seulement les gènes en commun. La valeur par défaut est TRUE.
-	`tables` : Une chaîne de caractère contenant les chemins vers les matrices de données normalisées à fusionner, séparés par des virgules. La valeur par défaut est "MicroArray_simulation.txt,RNA-seq_simulation.txt"

####    2.4.2	Plus de détails 

Un prétraitement est réalisé afin de rendre les données pertinentes et exploitables comme l’existence des doublons dans les noms de patients entre les tableaux. En effet, S'il y a des doublons, on ajoute un suffixe "_ti" au bout du nom.

Aussi, la vérification de la nature des données. Ils doivent être numériques, dans le cas contraire, ils seront transformés en données numériques. Puis une vérification des données manquantes est effectuée. Dans ce cas-là, la fonction s’arrête et renvoie un message d’erreur.

####    2.4.3	Sortie

La fonction renvoie une matrice de données fusionnées traitées avec la méthode de standardisation choisie.

####      2.4.4	Execution avec Rscript
```
Rscript rnaseq_microarray_fusion.r --standardisation zscore --tables MicroArray_simulation.txt, RNASeq_simulation.txt
```

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
R packages :
- optparse : https://CRAN.R-project.org/package=optparse
- edgeR : https://www.bioconductor.org/packages/release/bioc/html/edgeR.html
- DESeq2 : https://www.bioconductor.org/packages/release/bioc/html/DESeq2.html
- preprocessCore : https://www.bioconductor.org/packages/release/bioc/html/preprocessCore.html

