# Uses resnik and lin similiarities, rescaled between 0 and 1 and averaged.
# manual assignment MeSH expected to be formatted a la Nelson (lower case)
# can do a square matrix or (if don't set trait_terms) or set terms
# to be in rows and trait_terms to be in columns.

construct_MeSH_similarity <- function(terms, mesh_names, mesh_tree, manual_assignments, IC=NULL, trait_terms = terms) {
  mesh_df <- left_join(mesh_names %>% filter(Preferred=="Y"), mesh_tree)
  all_ontology_terms <- terms[terms %in% mesh_df$Name]
  all_ontology_traits <- trait_terms[trait_terms %in% mesh_df$Name]
  # Term similarity for the mesh
  average_similarity <- mesh_average_similarity(all_ontology_terms, mesh_df, IC, col_terms = all_ontology_traits)
  # Incorporate Nelson's manual similarities
  average_similarity <- add_manual_assignments(average_similarity, manual_assignments)
  # Get the terms that are not MeSH and assign zero similarity to all other terms
  non_mesh_terms <- terms[!terms %in% mesh_df$Name]
  non_mesh_traits <- trait_terms[!trait_terms %in% mesh_df$Name]
  non_mesh_rows <- matrix(data=0, nrow=length(non_mesh_terms), ncol=ncol(average_similarity) + length(non_mesh_traits))
  rownames(non_mesh_rows) <- non_mesh_terms
  non_mesh_cols <- matrix(data=0, nrow=nrow(average_similarity), ncol=length(non_mesh_traits))
  colnames(non_mesh_cols) <- non_mesh_traits
  average_similarity <- rbind(cbind(average_similarity, non_mesh_cols), non_mesh_rows)
  # If any of the elements are self-similarities, set to 1
  average_similarity <- set_self_similarity(average_similarity)
  return(average_similarity)
}

# Adding manual assignments (formated as in Nelson supplementary table with columns MSH1, MSH2, ManSim) to term similarity matrix
add_manual_assignments <- function(term_similarity_mat, manual_assignments) {
  all_terms <- unique(c(rownames(term_similarity_mat), colnames(term_similarity_mat)))
  manual_assignments <- mutate(manual_assignments, MSH1anum=to_alphanum(MSH1), MSH2anum=to_alphanum(MSH2))
  manual_assignments_red <- manual_assignments %>% 
    filter(MSH1anum %in% to_alphanum(all_terms), MSH2anum %in% to_alphanum(all_terms)) %>%
    left_join(data.frame(MSH10=all_terms, MSH1anum=to_alphanum(all_terms), stringsAsFactors = FALSE)) %>% 
    left_join(data.frame(MSH20=all_terms, MSH2anum=to_alphanum(all_terms), stringsAsFactors=FALSE)) %>% 
    filter(!is.na(MSH10), !is.na(MSH20))
  manual_assignments_red1 <- filter(manual_assignments_red, MSH10 %in% rownames(term_similarity_mat), MSH20 %in% colnames(term_similarity_mat))
  manual_assignments_red2 <- filter(manual_assignments_red, MSH20 %in% rownames(term_similarity_mat), MSH10 %in% colnames(term_similarity_mat))
  term_similarity_mat[cbind(manual_assignments_red1$MSH10, manual_assignments_red1$MSH20)] <- manual_assignments_red1$ManSim
  # Symmetry!
  term_similarity_mat[cbind(manual_assignments_red2$MSH20, manual_assignments_red2$MSH10)] <- manual_assignments_red2$ManSim  
  return(term_similarity_mat)
}

# compute average of lin and scaled resnik similarity of terms
mesh_average_similarity <- function(terms, mesh_df, IC=NULL, col_terms=terms) {
  if (!all(terms %in% mesh_df$Name)) {
    stop("Expect all terms to be MeSH heading")
  }
  # Using the ontology packages to compute similarities
  oi <- mesh_to_ontology_index(mesh_df)
  # descendants IC unless otherwise specified
  if (is.null(IC)) {
    IC <- descendants_IC(oi)
  }
  # Lin similarity and unnormalized Resnik similarity
  lin_similarity <- get_term_sim_mat(oi, information_content=IC, method = "lin", row_terms = terms, col_terms = col_terms)
  resnik_similarity <- get_term_sim_mat(oi, information_content=IC, method = "resnik", row_terms = terms, col_terms = col_terms)
  # Standardizing Resnik similarity
  resnik_similaritys <- resnik_similarity/max(resnik_similarity)
  # Set self similarities to 1.
  resnik_similaritys <- set_self_similarity(resnik_similaritys)
  # Average similarities a la Nelson
  average_similarity <- (resnik_similaritys + lin_similarity)/2
  return(average_similarity)
}

# Function to set self-similarities (same row name and column name) to 1.
set_self_similarity <- function(similarity_matrix) {
  matching_col = match(rownames(similarity_matrix), colnames(similarity_matrix))
  inds <- cbind((1:nrow(similarity_matrix))[!is.na(matching_col)], matching_col[!is.na(matching_col)])
  similarity_matrix[inds] <- 1
  return(similarity_matrix)
}

# Convert MeSH to ontologyIndex
mesh_to_ontology_index <- function(mesh_df) {
  mesh_headings_tree <- filter(mesh_df, Preferred=="Y", !is.na(TreeNumber))
  # Top level headings (C=disease, etc) are not included in this file but are part of the structure
  top_headings <- unique(substr(mesh_headings_tree$TreeNumber, 1, 1))
  mesh_headings_tree <- rbind(mesh_headings_tree, data.frame(UI=NA, Name=top_headings, Preferred="Y", TreeNumber=top_headings))
  # Parents of each node from tree number
  mesh_headings_tree <- mutate(mesh_headings_tree, 
                               Parent=sapply(strsplit(TreeNumber, "\\."), 
                                             function(x) ifelse(length(x) > 1, paste(x[-length(x)], collapse = "."), 
                                                                ifelse(nchar(x)==1, "MeSH", substr(x,1,1)))))
  # Name for each parent
  mesh_headings_tree <- left_join(mesh_headings_tree, mesh_headings_tree %>% dplyr::select(ParentName=Name, TreeNumber),
                                  by=c("Parent"="TreeNumber"))
  mesh_headings_tree <- replace_na(mesh_headings_tree, list(ParentName="MeSH"))
  # List of parents of each node
  parents_list <- split(mesh_headings_tree$ParentName, mesh_headings_tree$Name)
  parents_list <- c(parents_list, list(MeSH=c()))
  # Ontology index
  oi = ontology_index(parents_list)
  return(oi)
}

# Alphanumeric lowercase only, used for matching terms
to_alphanum <- function(x) {
  tolower(gsub("[^[:alnum:]]", "", x))
}