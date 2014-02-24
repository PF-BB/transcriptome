Transcriptome
=============

Analyse de données transcriptomiques.

Contient : 
  - le package d'analyse SIMA (SImple Microarray Analysis).
  - une ébauche d'application shiny pour visualiser le profil d'expression d'un gène dans plusieurs jeux de données à la fois.

### Nomenclature des fichiers
  - ``01_`` :arrow_right: normalisation.
  - ``02_`` :arrow_right: analyse différentielle.
  - ``03_`` :arrow_right: plots.
  - ``04_`` :arrow_right: annotation.
  - ``05_`` :arrow_right: écriture des tables.

### Pour tester SIMA sans re-générer le Package

```S
rm(list=ls())
setwd("~/git/transcriptome")

require(devtools)
require(roxygen2)

# load SIMA (dev) 
load_all("SIMA")

# generate documentation
document("SIMA")
```
