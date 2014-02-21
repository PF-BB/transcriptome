##########################################################
## Analyse de donnees Microarray avec la librairie SIMA ##
##########################################################

## 0. Preambule

# Nettoyer la memoire
rm(list=ls())

# Charger la librairie SIM
library(SIMA)

# Declarer les chemins
# --> ou se trouvent les donnees microarray + les phenotypes d'interet
DATA_PATH <- "~/data" # 
# --> ou ecrire les resultats
OUT_PATH <- "~/resultats"

## 1. Normalisation
eset.norm <- normalisation(in="")

## 2. Analyse differentielle
contr <- ""
res <- analyse_diff(eset,contr)

## 3. Representations graphiques
acp(eset)
cluster(eset)
volcanoplot(res)

## 4. Annotation
eset.annot <- annotation(eset)

## 5. Ecriture des tables de resultats
ecrire(res)