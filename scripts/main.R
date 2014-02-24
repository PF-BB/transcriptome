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
DATA_PATH <- "~/data" # ou bien ? ILLUMINA_FILE_PATH <- "~/data/raw_data_illumina.txt"
# --> ou ecrire les resultats
OUT_PATH <- "~/resultats"

preambule(DATA_PATH, OUT_PATH)

## 1. Normalisation
eset.norm <- normalisation(in=DATA_PATH)

## 2. Analyse differentielle
y <- read.table( file.path(DATA_PATH,"phenotype.txt") , header = TRUE)[,2]
contr <- declare_contrastes(y)
res <- analyse_diff(eset, y, contr)

## 3. Representations graphiques
acp(eset)
cluster(eset)
volcanoplot(res)

## 4. Annotation
eset.annot <- annotation(eset)

## 5. Ecriture des tables de resultats
ecrire(res)