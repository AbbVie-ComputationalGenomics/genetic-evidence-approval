# genetic-evidence-approval
Code and supplementary to reproduce main text figures on new estimates of the effect of genetic evidence on drug approval.
## Quickstart guide
This guide gives instructions for compiling a pdf document reproducing main text figures from King et al. 2018 using supplementary data sources and for installing the associated shiny app outputing predictions.
### Cloning repository
`git clone https://pig.abbvienet.com/kingea/genetic-evidence-approval/kingea/genetic-evidence-approval.git`

Note: For public git repo, to appear,

`git clone https://github.com/AbbVie-ComputationalGenomics/genetic-evidence-approval.git`

From an R session
##### Figure reproduction dependencies
```
install.packages("knitr")
install.packages("dplyr")
install.packages("ggplot2")
install.packages("cowplot")
install.packages("gdata")
install.packages("epitools")
install.packages("biomaRt")
```
##### Shiny app dependencies
```
install.packages("shiny")
install.packages("DT")
install.packages("dplyr")
```
##### Additional dependencies for running stan models
```
install.packages("rstan")
```
##### Additional dependencies for ontology
```
install.packages("ontologySimilarity")
install.packages("ontologyIndex")
```
### Downloading data and results files
This is how to download model output if not running Stan.
```
cd path/to/genetic-evidence-approval
mkdir data
mkdir results
cd results
[fetch ShinyAppPrecomputed.rds from paper supplementary materials]
```
 Or download directly from paper and move to `path/to/genetic-evidence-approval/results`
#### File only data files
##### Nelson et al supplementary tables
From command line
```
cd path/to/genetic-evidence-approval/data
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S12.txt
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S13.txt
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S14.txt
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S9.xlsx
 ```
 Or download directly from paper and move to `path/to/genetic-evidence-approval/data`
 ##### King et al supplementary tables
 ```
cd path/to/genetic-evidence-approval/data
[fetch gene_trait_assoc.tsv, Standardized_Nelson_Associations.tsv, target_indication_nmsh.tsv, target_indication.tsv, Target_Properties.tsv from paper supplementary materials]
```
 Or download directly from paper and move to `path/to/genetic-evidence-approval/data`
### Reproducing main text figures
From command line
```
Rscript -e "library(knitr); knit('doc/AssociationBetweenGeneticEvidenceAndSuccess.Rnw')"
pdflatex doc/AssociationBetweenGeneticEvidenceAndSuccess.tex
```
From RStudio, open `AssociationBetweenGeneticEvidenceAndSuccess.Rnw` and click compile pdf.  Ensure Sweave options are set to weave files using knitr.
[TODO: convert this to a markdown document to avoid complications]
### Rerunning model fit
```
Rscript -e "library(knitr); knit('doc/StanModelFits.Rnw')"
```
From RStudio, open `StanModelFits.Rnw` and click compile pdf.  Ensure Sweave options are set to weave files using knitr.  This is likely to take hours to run.
### Running shiny app
From RStudio, open `PredictApproval/app.R` and click run app button.
[TODO how to from Rsession]
[TODO additional shiny app instructions]

