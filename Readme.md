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
### 	microarray_simul.r
...

## Deployment and usage in Galaxy workflows
...

## Software dependencies
