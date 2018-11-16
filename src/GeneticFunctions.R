get_xMHC_genes <- function(xMHC_start_ensembl='ENSG00000164508', xMHC_end_ensembl='ENSG00000237649', max_try=10) {
  require(biomaRt)
  ntry <- 0
  repeat {
    ntry <- ntry + 1
    mart <- try(useMart("ENSEMBL_MART_ENSEMBL"))
    if (!inherits(mart, "try-error")) {
      human_genes <- try(useDataset("hsapiens_gene_ensembl", mart))
      if (!inherits(human_genes, "try-error")) {
        ensembl_lookup <- try(getBM(attributes=c('ensembl_gene_id','external_gene_name','entrezgene','chromosome_name', 'start_position', 'end_position'), 
                                    filters='chromosome_name', values = 6, mart = human_genes))
        if (!inherits(ensembl_lookup, "try-error")) break
      }
      if (ntry==max_try) break
    }
  }
  xMHC_start <- filter(ensembl_lookup, ensembl_gene_id==xMHC_start_ensembl) %>% pull(start_position)
  xMHC_end <- filter(ensembl_lookup, ensembl_gene_id==xMHC_end_ensembl) %>% pull(end_position)
  xMHC_ensembl <- filter(ensembl_lookup, !is.na(ensembl_gene_id), 
                                   start_position >= xMHC_start, end_position <= xMHC_end) %>% pull(ensembl_gene_id)
  return(xMHC_ensembl)
}
