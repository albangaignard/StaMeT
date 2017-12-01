# Multi-transcripts toolbox

## Description
...

## Motivations
...


I.	La simulation 

Afin de pouvoir utiliser au mieux les outils d’analyses, quelques légères modifications ont été apportées aux fonctions de simulation prises en référence et ce, pour chacune des technologies. Ces modifications concernaient principalement l’organisation des matrices de données obtenues en sortie.

1.	microarray_simul.r

Pour simuler les données microarrays, nous avons choisi l’algorithme ‘’Madsim’  (Réf), implémenté en R, nous avons ajouté deux autres paramètres permettant de contrôler la variation du vecteur des fold-changes (FC).

1.1.1	Description

Cette fonction permet de simuler des données microarrays à partir d’un modèle prédéfini. Les données simulées ont un comportement similaire aux données microarrays produites par la plateforme « Affimetrix ».
Pour cela, l’utilisateur doit fournir un ensemble des paramètres, ou utiliser ceux disponibles par défaut.

1.1.2	Usage

" microarray_simul.r" --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 –m1 1.4 --m2 0.8

1.1.3	Arguments

    •	"-gn" ou  "--gene_number": Un nombre entier naturel indiquant  le nombre de gènes dans les données simulées. La valeur par défaut est : --gene_number=10,000
    •	"-sn1" ou  "--samples_n1" : Un nombre entier naturel indiquant le nombre d’échantillons  du phénotype 1 (condition 1). La valeur par défaut est : --samples_n1=75
    •	"-sn2" ou  "--samples_n2" : Un nombre entier naturel indiquant le nombre d’échantillons  du phénotype 1 (condition 1). La valeur par défaut est : --samples_n2=75
    •   "-diff "ou  " --diff_genes_ratio" : Un nombre décimal indiquant le the pourcentage des gènes différentiellement exprimés. Sa valeur par défaut est : --diff_genes_ratio=0.1
    •	"-up", ou "--up_ratio" : Un nombre décimal indiquant le pourcentage de gènes surexprimés. Sa valeur par défaut est : --up_ratio=0.5
    •	"-m1" ou "--m1": Un nombre décimal  correspondant à la différence moyenne entre la moyenne totale et la moyenne  des gènes différentiellement exprimé avec  des valeurs élevées. Sa valeur par défaut est : --m1=1.4
    •	"-m2" ou  "--m2": Un nombre décimal  correspondant à la différence moyenne entre la moyenne totale et la moyenne  des gènes différentiellement exprimé avec  des valeurs peu elevées. Sa valeur par défaut est : --m1=0.8
    •	"-s" ou "--seed": Un entier utilisé pour générer un nombre aléatoire par l'ordinateur



1.1.4	Plus de détails :

-- Si l’utilisateur fournit un nombre décimal au lieu d’un nombre  entier  pour les trois premiers paramètres, la valeur sera arrondie.

-- La fonction ne sera pas exécutée et retournera un message d’erreur dans les cas suivant :

    •	Si un des paramètres numérique ne l’est pas

    •	Si un des paramètres numériques est négatif 

    •	Si le nombre de gènes à simuler est nul

    •	Si le nombre d’échantillons est nul 

    •	Si les deux paramètres de proportions ne sont pas compris entre 0 et 1



1.1.5	Sortie 

La fonction renvoie une matrice de données avec respectivement  le nombre de lignes et de colonnes spécifié par les paramètres d'entrée "--gene_number et "--samples_n1 + "--samples_n2. 
Les données sont supposées  etre semblables aux données microarrays produites par la plateforme « Affimetrix » log2 intensité.

1.1.6	Exemple 

RScript " microarray_simul.r" --gene_number 1000 --samples_n1 20 --samples_n2 20 --up_ratio 0.5 --diff_genes_ratio 0.1 –m1 1.4 --m2 0.8
