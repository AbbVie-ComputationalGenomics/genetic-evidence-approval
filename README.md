# genetic-evidence-approval
Code and supplementary to reproduce main text figures on new estimates of the effect of genetic evidence on drug approval.
## Quickstart guide
This guide gives instructions for compiling a pdf document reproducing main text figures from King et al. 2018 using supplementary data sources and for installing the associated shiny app outputing predictions.
### Cloning repository
`git clone https://pig.abbvienet.com/kingea/genetic-evidence-approval.git`

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

#### File only data files
##### Nelson et al supplementary tables
From command line
```
cd genetic-evidence-approval
mkdir data
cd data

wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S12.txt
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S13.txt
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S14.txt
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S9.xlsx
 ```
 Or download directly from paper and move to `path/to/genetic-evidence-approval/data`
 ##### King et al supplementary tables
 ```
[fetch gene_trait_assoc.tsv, Standardized_Nelson_Associations.tsv, target_indication_nmsh.tsv, target_indication.tsv, Target_Properties.tsv from paper supplementary materials]
```
 Or download directly from paper and move to `genetic-evidence-approval/data`
### Reproducing main text figures
From command line
```
Rscript -e "library(knitr); knit('doc/AssociationBetweenGeneticEvidenceAndSuccess.Rmd')"
```
From RStudio, open `AssociationBetweenGeneticEvidenceAndSuccess.Rmd` and click Knit.
### Rerunning model fit
From command line
```
Rscript -e "library(knitr); knit('doc/StanModelFits.Rmd')"
```
From RStudio, open `StanModelFits.Rmd` and click Knit.
### Running shiny app
From RStudio, open `PredictApproval/app.R` and click run app button.
[TODO how to from Rsession]
[TODO additional shiny app instructions]

