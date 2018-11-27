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
install.packages("markdown")
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
##### Nelson et al supplementary tables and King et al supplementary tables
From command line
```
cd genetic-evidence-approval
mkdir data
cd data

wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S12.txt
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S13.txt
wget https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S14.txt

[fetch gene_trait_assoc.tsv, Standardized_Nelson_Associations.tsv, target_indication_nmsh.tsv, target_indication.tsv, Target_Properties.tsv, top_mesh.tsv, indication_trait_similarity.tsv from paper supplementary materials.  url not yet available]

cd ../
```
 Or download directly from papers and move to `genetic-evidence-approval/data`
### Reproducing main text figures
From R session in `genetic_evidence_approval` directory
```
library(knitr)
library(markdown)
knit('doc/AssociationBetweenGeneticEvidenceAndSuccess.Rmd', 'doc/AssociationBetweenGeneticEvidenceAndSuccess.md')
markdownToHTML('doc/AssociationBetweenGeneticEvidenceAndSuccess.md', 'doc/AssociationBetweenGeneticEvidenceAndSuccess.html')
browseURL(paste('file://', file.path('doc/AssociationBetweenGeneticEvidenceAndSuccess.html'), sep=''))
```
From RStudio, open `AssociationBetweenGeneticEvidenceAndSuccess.Rmd` and click Knit.

### Rerunning model fit
From R session in `genetic_evidence_approval` directory
```
library(knitr)
library(markdown)
knit('doc/StanModelFits.Rmd', 'doc/StanModelFits.md')
markdownToHTML('doc/StanModelFits.md', 'doc/StanModelFits.html')
browseURL(paste('file://', file.path('doc/StanModelFits.html'), sep=''))
```
From RStudio, open `StanModelFits.Rmd` and click Knit.
Note this will likely take several hours to run.
### Running shiny app
From R session in `genetic_evidence_approval` directory
```
library(rshiny)
runApp('PredictApproval/app.R')
```
From RStudio, open `genetic_evidence_approval/PredictApproval/app.R` and click run app button.

## Notes on shiny app
The shiny app displays the estimated success probability of gene target-indication pairs, and the odds ratio of success given genetic evidence.  
### Using the app
##### Tabs
The app can display genetic evidence by target for a fixed indication or genetic evidence by indication for a fixed gene target.  These are two separate tabs.
##### Available models
GWAS: GWAS genetic evidence alone
OMIM: OMIM genetic evidence alone
GWAS and OMIM: GWAS and OMIM genetic evidence in the same model.  
##### Available gene target-indication pairs
Gene target indication pairs were created from MeSH terms mapping to Pharmaprojects indications and genetically associated non xMHC, protein coding genes (filtering criteria used in the paper).  Target-indication pairs are only included if the target is genetically linked to a trait with similarity at least 0.5.  
##### Output interpretation
Success probabilities reflect target and indication level properties in addition to genetic evidence, and are best interpreted on a relative scale due to unknown development times for new programs (see next section for further details).  Odds ratios purely reflect the contribution of genetic evidence to success, and by default we sort by the odds ratio.  Odds ratios above 1 mean the gene target-indication pair is more likely to succeed than if there were no genetic association and odds ratios below 1 mean it is less likely to succeed than with no genetic assocation.  
### Interpretation caveats
##### Prior plausibility
Note that fitted probabilities are produced from Pharmaprojects gene target-indication pairs, which had sufficient prior plausibility to enter early stage (Preclinical or Phase I) development.  Therefore predicted probabilities may be inflated for gene target-indication pairs with low prior plausibility (for example, due to negative early discovery results or problems with druggability) as these were not part of the training data.  The option to restrict to known targets when searching for targets by indication is intended to reduce this problem.  
##### Developpment Time
Fitted probabilities use a fixed value for time target has been under development for all targets.  This is because we did not want estimates to reflect how much work has already been done on the target for the purposes of comparing genetic evidence, and because the development time for new targets is not known.  Because of this probabilities are better interpreted on a relative than absolute scale.

