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
install.packages("tidyr")
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
From R session in genetic_evidence_approval/doc directory
```
download.file('https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S12.txt', '../data/ng.3314-S12.txt')
download.file('https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S13.txt', '../data/ng.3314-S13.txt')
download.file('https://images.nature.com/full/nature-assets/ng/journal/v47/n8/extref/ng.3314-S14.txt', '../data/ng.3314-S14.txt')
```
 Or download directly from papers and move to `genetic-evidence-approval/data`
### Reproducing main text figures
From R session in `genetic_evidence_approval/doc` directory
```
library(knitr)
library(markdown)
knit('AssociationBetweenGeneticEvidenceAndSuccess.Rmd', 'AssociationBetweenGeneticEvidenceAndSuccess.md')
markdownToHTML('AssociationBetweenGeneticEvidenceAndSuccess.md', 'AssociationBetweenGeneticEvidenceAndSuccess.html')
browseURL(paste('file://', file.path(getwd(),'AssociationBetweenGeneticEvidenceAndSuccess.html'), sep='')) 
```
From RStudio, open `AssociationBetweenGeneticEvidenceAndSuccess.Rmd` and click Knit.  The file will take several minutes to run.

### Rerunning model fit
From R session in `genetic_evidence_approval/doc` directory
```
library(knitr)
library(markdown)
knit('StanModelFits.Rmd', 'StanModelFits.md')
markdownToHTML('StanModelFits.md', 'StanModelFits.html')
browseURL(paste('file://', file.path(getwd(),'StanModelFits.html'), sep='')) 
```
From RStudio, open `StanModelFits.Rmd` and click Knit.
Note this will likely take several hours to run.
### Running shiny app
From R session in `genetic_evidence_approval` directory
```
library(shiny)
runApp('PredictApproval/app.R')
```
From RStudio, open `genetic_evidence_approval/PredictApproval/app.R` and click run app button.

## Notes on shiny app
The shiny app displays the estimated success probability of gene target-indication pairs, and the odds ratio of success given genetic evidence.  
### Using the app
##### Tabs
The app can display genetic evidence by target for a fixed indication or genetic evidence by indication for a fixed gene target.  These are two separate tabs.
##### Available models
* GWAS: GWAS genetic evidence alone
* OMIM: OMIM genetic evidence alone
* GWAS and OMIM: GWAS and OMIM genetic evidence in the same model.  
##### Available gene target-indication pairs
Gene target indication pairs were created from MeSH terms mapping to Pharmaprojects indications and genetically associated non xMHC, protein coding genes (filtering criteria used in the paper).  Target-indication pairs are only included if the target is genetically linked to a trait with similarity at least 0.5 to some indication.  
##### Output interpretation
Success probabilities reflect target and indication level properties in addition to genetic evidence, and are best interpreted on a relative scale due to unknown development times for new programs (see next section for further details).  Odds ratios purely reflect the contribution of genetic evidence to success, and by default we sort by the odds ratio.  Odds ratios above 1 mean the gene target-indication pair is more likely to succeed than if there were no genetic association and odds ratios below 1 mean it is less likely to succeed than with no genetic assocation.  
### Interpretation caveats
##### Prior plausibility
Note that fitted probabilities are produced from Pharmaprojects gene target-indication pairs, which had sufficient prior plausibility to enter early stage (Preclinical or Phase I) development.  Therefore predicted probabilities may be inflated for gene target-indication pairs with low prior plausibility (for example, due to negative early discovery results or problems with druggability) as these were not part of the training data.  The option to restrict to known targets when searching for targets by indication is intended to reduce this problem.  
##### Development Time
Fitted probabilities use a fixed value for time target has been under development for all targets.  This is because we did not want estimates to reflect how much work has already been done on the target for the purposes of comparing genetic evidence, and because the development time for new targets is not known.  Because of this probabilities are better interpreted on a relative than absolute scale.

## Supplementary data files
###### target_indication.tsv
Data file with one row per MSH-ensembl_id pair derived from Informa Pharmaprojects.
* lApprovedUS.EU = is there a US/EU approved drug with this target approved for this indication?  
* Phase.Latest = Inferred latest historical development phase from Pharmaprojects (see methods and supplement for construction).  
* First Added = earliest added date for any drug with this target and indication.  
* Inactive = Is any drug with this target and indication assigned an inactive status in Pharmaprojects, such as No Development Reported?  * symbol = HUGO/HGNC symbol corresponding to the target ensembl id.
##### target_indication_nmsh.tsv
Same as above, but the MeSH term mappings have been harmonized with Nelson et al. 2015.
##### gene_trait_assoc.tsv
GWAS catalog and OMIM gene-trait links used in this analysis.  Some columns are only applicable to GWAS Catalog or OMIM but not both, and some GWAS Catalog columns are only applicable for some types of evidence.  Only one LD SNP is retained per top SNP-gene-trait triplet, and which is retained is determined by the highest score (unless there is a moderate/high deleterious variant with higher LD than the top scoring: we need to retain these for the analysis of high confidence SNP-Gene links).  Key columns are: 
* SNP_A: the GWAS catalog top SNP
* SNP_B the LD SNP providing the information in the row
* R2 = LD between SNP_A and SNP_B
* ensembl_id = ensembl_id of associated gene
* eQTL = whether or not SNP_B is a GTeX eQTL for ensembl_id (p < 10e-6 in some tissue)
* eQTL_pval = lowest p-value of any tissue in GTeX for SNP_B and ensembl_id
* DHS = is there an established DHS link (see methods)
* Distance = distance between gene and SNP_B if less than or equal to 5000 bp
* Deleterious = SnpEff assessed deleteriousness of SNP_B on gene, missense = is SNP_B a missense variant in gene
* distance = is SNP_B within 5000 bp of gene
* RDB = RegulomeDB score for SNP_B if 4 or less
* MAPPED_TRAIT = GWAS Catalog MAPPED_TRAIT for association
* TotalScore = score for association according to Nelson et al. method
* MSH = mapped MeSH term, Rank = rank of gene for that assocation (relative to other associated genes)
* NStudy = number of GWAS catalog studies finding SNP_A linked to MAPPED_TRAIT after filtering on p < 1e-8
* Source = GWAS:A for GWAS Catalog, OMIM for OMIM
* Phenotype = OMIM phenotype
* UI = MeSH unique identifier
* MSH = mapped MeSH heading
* first_added = earliest date a link between the MeSH heading and ensembl_id appeared in the source (this may differ from the date the association appeared because more than one association may have the same mapped MeSH term and ensembl_id)
* symbol = HUGO/HGNC symbol for gene
* xMHCGene = is the gene considered part of xMHC?
##### Standardized_Nelson_Associations.tsv
Nelson et al. supplementary table 2, standardized to use 2017 MeSH and ensembl gene ids for comparison with new data.
##### indication_trait_similarity.tsv
A matrix with MeSH terms mapped to Pharmaprojects indications as rows and MeSH terms mapped to genetically associated traits as columns with the entries being the similarity between the terms.
##### Target_properties.tsv
Data file with one row per ensembl_id giving 
* GO annotations for each target (1 = has annotation, 0 = does not).  GO annotations were a subset of all possible GO terms selected to map to Pharmaprojects target families.
* RVIS = residual variation intolerance score from [Petrovski 2013](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1003709) supplement
* Time = time drug target has been under development, in days since Jan 1, 1970.
##### top_mesh.tsv
For each MeSH heading, gives all top-level MeSH (e.g. Neoplasms) under which the heading appears.

