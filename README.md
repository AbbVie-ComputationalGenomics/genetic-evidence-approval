# genetic-evidence-approval
Code and supplementary to reproduce main text figures on new estimates of the effect of genetic evidence on drug approval.
## Quickstart guide
This guide gives instructions for ...
### Cloning repository
### Installing dependencies
##### Common dependencies
##### Shiny app only dependencies
##### File only dependencies
### Downloading data and results files
##### Common data files
[we do not have these files hosted anywhere and if on github they will need to be under git control.  Code to reproduce MeSH files???  Or just make them a supplementary file even though they have nothing to do with the paper?  Similarity, I don't want to make HGNC protein coding genes a supplementary material...I can provide my python script to parse mesh, they can download 2017 mesh from server and run that.  And I think they can download their own copy of HGNC coding]
```
cd path/to/genetic-evidence-approval
mkdir data
mkdir results
cd results
[fetch ShinyAppPrecomputed.rds from paper supplementary materials]
```
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
 ##### King et al supplementary tables
 ```
cd path/to/genetic-evidence-approval.data
[fetch gene_trait_assoc.tsv, Standardized_Nelson_Associations.tsv, target_indication_nmsh.tsv, target_indication.tsv, Target_Properties.tsv from paper supplementary materials]
```
### Reproducing main text figures
### Running shiny app
### Rerunning model fit
