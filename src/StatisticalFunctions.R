##################################################################
##### Nelson et al. Table 1 Replication ##########################
##################################################################

# A function that takes a drug table and annotates it with
# whether or not each target appears in OMIM or GWAS gene list
# and the closest trait similarity for supporting genetic association
# of each type.  
# Supporting trait similarity set to zero if target not associated.
# Default filter list is set to any genetic association (any), omim only, and gwas only, but
# this can be changed to allow annotation with different types of genetic evidence for different
# filtering functions (e.g. filter on missense variants)
# Should be named to label columns

annotate_target_indication_table_with_genetic_evidence <- function(target_indication_table, association_table, 
                                                                   MSH_similarity, gene_col_name, 
                                                                   filter_list=list(gwas=gwas_filter, omim=OMIM_filter, any=function(x) {x}),
                                                                   original_trait_cols = c("MAPPED_TRAIT", "Phenotype", NA)) {
#  gene_col2 = sym(gene_col_name)
  for (j in seq_along(filter_list)) {
    filtered_table <- filter_list[[j]](association_table)
    filter_name <- names(filter_list)[j]
    target_indication_table[[gene_intercept_colname(filter_name)]] <- target_indication_table[[gene_col_name]] %in% filtered_table[[gene_col_name]]
    target_indication_table[[similarity_base(filter_name)]] <- rep(0, nrow(target_indication_table))
    if (!is.na(original_trait_cols[j])) {
      target_indication_table[[trait_colname(filter_name)]] <- rep(NA, nrow(target_indication_table))
    }
    for (i in which(target_indication_table[[gene_intercept_colname(filter_name)]])) {
      target_gene <- target_indication_table[[gene_col_name]][i]
#      supporting_msh_rows <- filter(filtered_table, !!gene_col2 == target_gene)
      supporting_msh_rows <- filtered_table[filtered_table[[gene_col_name]]==target_gene,]
      supporting_msh <- supporting_msh_rows$MSH
      supporting_similarity <- MSH_similarity[target_indication_table$MSH[i],supporting_msh]
      target_indication_table[[similarity_base(filter_name)]][i] <- max(supporting_similarity)
      if (!is.na(original_trait_cols[j])) {
        target_indication_table[[trait_colname(filter_name)]][i] <- supporting_msh_rows[[original_trait_cols[j]]][which.max(supporting_similarity)]
      }
    }
  }
  return(target_indication_table)
}

gene_intercept_colname <- function(assoc_class) {
  paste0(assoc_class, "_gene") 
}

similarity_base <- function(assoc_class) {
  paste0("best_mesh_", assoc_class)
}

trait_colname <- function(assoc_class) {
  paste0("best_trait_", assoc_class)
}

# Function to replicate Table1 of Nelson et al 2015 from drug, association, and similarity data as well as cutoff parameters.
# MSH_to_include is an optional argument giving the set of MeSH terms to include in the analysis.  Otherwise those with at least
# ngene_cutoff associations are included.

replicate_table1 <- function(target_indication_table, association_table, MSH_similarity, 
                             similarity_cutoff, ngene_cutoff, gene_col_name, MSH_to_include=NULL, 
                             source_names=c("gwas"="GWASdb", "omim"="OMIM", "any"="GWASdb & OMIM")) {
  if (is.null(MSH_to_include)) {
    MSH_to_include <- get_well_studied_MSH(target_indication_table, association_table, MSH_similarity, similarity_cutoff, ngene_cutoff, gene_col_name)
  }
  filtered_target_indication_table <- filter(target_indication_table, MSH %in% MSH_to_include)
  annotated_target_indication_table <- 
    annotate_target_indication_table_with_genetic_evidence(filtered_target_indication_table, association_table, MSH_similarity, gene_col_name, original_trait_cols = rep(NA,3))
  annotated_target_indication_table <- filter_and_label_by_phase(annotated_target_indication_table)
  # Trial pairs and sources for which risk ratios to be constructed
  trial_pairs <- cbind(c("PhaseI", "PhaseII", "PhaseIII", "PhaseI", "PhaseI"),
                       c("PhaseII", "PhaseIII", "Approved", "PhaseIII", "Approved"))
  sources <- c("gwas", "omim", "any")
  df <- vector("list", length(sources))
  #risk ratio tests.
  for (i in 1:length(sources)) {
    source_df <- vector("list", length(trial_pairs))
    for (j in 1:nrow(trial_pairs)) {
      source_df[[j]] <- phase_transition_relative_risk_df(annotated_target_indication_table, trial_pairs[j,1], trial_pairs[j,2], sources[i], similarity_cutoff, ngene_cutoff)
    }
    df[[i]] <- bind_rows(source_df)
  }
  df <- bind_rows(df)
  df$FirstPhase <- gsub("PhaseI", "Phase I", df$FirstPhase)
  df$SecondPhase <- gsub("PhaseI", "Phase I", df$SecondPhase)
  df <- mutate(df, Source=source_names[Source])
  return(df)
}

# Filter to retain only new trial phases and label with whether each target-indication pair has attained at least phase x
filter_and_label_by_phase <- function(target_indication_table) {
  interesting_trial_phases <- c("Phase I Clinical Trial", "Phase II Clinical Trial", "Phase III Clinical Trial")
  target_indication_table <- filter(target_indication_table, Phase.Latest %in% interesting_trial_phases | lApprovedUS.EU)
  target_indication_table <- mutate(target_indication_table, PhaseI=TRUE, 
                                    PhaseII=Phase.Latest %in% c("Phase II Clinical Trial", "Phase III Clinical Trial") | lApprovedUS.EU,
                                    PhaseIII=Phase.Latest=="Phase III Clinical Trial" | lApprovedUS.EU,
                                    Approved=lApprovedUS.EU)
  return(target_indication_table)
}

# one data frame row for a phase transition
phase_transition_relative_risk_df <- function(annotated_target_indication_table, first_phase, second_phase, source, similarity_cutoff, ngene_cutoff) {
  rrout <- progression_test(annotated_target_indication_table, first_phase, second_phase, source, similarity_cutoff)  
  df_entry <- data_frame_entry_from_risk_ratioboot(rrout, first_phase, second_phase, source, similarity_cutoff, ngene_cutoff)
  df_entry$N1 <- sum(annotated_target_indication_table[[first_phase]])
  df_entry$N2 <- sum(annotated_target_indication_table[[second_phase]] & annotated_target_indication_table[[first_phase]])
  return(df_entry)
}

#first_phase and second_phase are names of logical columns in annotated_target_indication_table

progression_test <- function(annotated_target_indication_table, first_phase, second_phase, source, similarity_cutoff) {
  at_least_second <- annotated_target_indication_table[[first_phase]] & annotated_target_indication_table[[second_phase]]
  at_first_stop_before_second <- annotated_target_indication_table[[first_phase]] &
    !at_least_second
  # TODO: switch to quo/enquo/!!
  has_association <- annotated_target_indication_table[[paste0("best_mesh_", source)]] >= similarity_cutoff
  nAN <- sum(has_association & at_first_stop_before_second)
  nAP <- sum(has_association & at_least_second)
  nNN <- sum(!has_association & at_first_stop_before_second)
  nNP <- sum(!has_association & at_least_second)
  twobytwo <- t(matrix(data=c(nNN, nNP, nAN, nAP), nrow=2, ncol=2))
  riskratio.boot(twobytwo)
}

# Function to construxt a row of my table 1 output from risk ratio book results.
data_frame_entry_from_risk_ratioboot <- function(rrboot_out, first_phase, second_phase, source, similarity_cutoff, ngene_cutoff, rrboot_out_row="Exposed2") {
  if (inherits(rrboot_out, "try-error")) {
    rrow <- c(estimate=NA, lower=NA, upper=NA)  
  } else {
    rrow <- rrboot_out$measure[rrboot_out_row,]
  }
  df <- data.frame(FirstPhase=first_phase, SecondPhase=second_phase, Source=source
                   , Estimate=rrow["estimate"], Lower=rrow["lower"], Upper=rrow["upper"],
                   SimilarityCutoff=similarity_cutoff, NGeneCutoff=ngene_cutoff, stringsAsFactors = FALSE)
  return(df)
}

# Filter association table to get OMIM only
# Excluding Omim is to improve agreement with main text of Nelson et al. 2015.
# All OMIM associations in this study are called OMIM.

OMIM_filter <- function(assoc_table) {
  require(dplyr)
  #  filter(assoc_table, Source %in% c("Omim", "OMIM"))  
  filter(assoc_table, Source=="OMIM")
}

# Filter association table to get any GWAS
# Including Omim produces results consistent with main text of Nelson et al. 2015.

gwas_filter <- function(assoc_table) {
  require(dplyr)
  #  filter(assoc_table, !(Source %in% c("Omim", "OMIM")))
  dplyr::filter(assoc_table, !Source=="OMIM")
}

# Get the MeSH terms with at least ngene_cutoff associations for similar traits.

get_well_studied_MSH <- function(target_indication_table, association_table, MSH_similarity, similarity_cutoff, ngene_cutoff, gene_col_name) {
  MSH_nassoc <- get_mesh_nassoc(target_indication_table, association_table, MSH_similarity, similarity_cutoff, gene_col_name)
  MSH_to_include <- filter(MSH_nassoc, NAssoc >= ngene_cutoff) %>% pull(MSH)
  return(MSH_to_include)
}

# Number of associations for each mesh term in target_indication_table for similar traits

get_mesh_nassoc <- function(target_indication_table, association_table, MSH_similarity, similarity_cutoff, gene_col_name) {
  umsh <- unique(target_indication_table$MSH)
  nassoc <- integer(length(umsh))
  for (i in 1:length(umsh)) {
    sim_MSH <- colnames(MSH_similarity)[MSH_similarity[umsh[i],] >= similarity_cutoff]
    association_rows <- filter(association_table, MSH %in% sim_MSH)
    nassoc[i] <- count_associations(association_rows)
  } 
  res <- data.frame(MSH=umsh, NAssoc=nassoc, stringsAsFactors = FALSE)
  return(res)
}

# Number of associations for each associated trait
# A distinct association is: distinct entry + MSH (OMIM or GWAS catalog) + distinct SNP id (for GWAS catalog)

count_associations <- function(association_rows) {
  return(length(unique(paste(association_rows$Link, association_rows$snp_id, association_rows$MSH))))
}

##################################################################
##### Progression Test ###########################################
##################################################################

# target_indication should be data frame with columns ensembl_id-MSH-Phase.Latest-lApprovedUS.EU as before with additional column NelsonStatus
generate_table1_progression <- function(target_indication_table, association_table, MSH_similarity, 
                                        similarity_cutoff, ngene_cutoff=5, gene_col_name="ensembl_id",
                                        source_names=c("gwas"="GWASdb", "omim"="OMIM", "any"="GWASdb & OMIM")) {
  well_studied <- get_well_studied_MSH(target_indication_table = target_indication_table, 
                                       association_table = association_table, 
                                       MSH_similarity = MSH_similarity, 
                                       similarity_cutoff = similarity_cutoff, ngene_cutoff = ngene_cutoff)
  candidates_to_advance <- filter(target_indication_table, 
                                  grepl("Trial", NelsonStatus), grepl("Clinical|Approved|egistr", Phase.Latest) | lApprovedUS.EU, 
                                  MSH %in% well_studied_clinical) %>% 
    mutate(Highest=pairwise_top_status(Phase.Latest, NelsonStatus), Valid=Highest==Phase.Latest) %>% 
    filter(Valid) %>% 
    mutate(Progressed=!Highest==NelsonStatus)
  annotated_candidates_to_advance <-
    annotate_target_indication_table_with_genetic_evidence(target_indication_table = candidates_to_advance,
                                                           association_table = association_table, MSH_similarity = MSH_similarity, 
                                                           gene_col_name = gene_col_name, original_trait_cols = rep(NA, 3))
  annotated_candidates_to_advance <- annotated_candidates_to_advance %>%
    mutate(PhaseIN=NelsonStatus=="Phase I Clinical Trial", 
           PhaseIIN=NelsonStatus=="Phase II Clinical Trial", 
           PhaseIIIN=NelsonStatus=="Phase III Clinical Trial",
           PhaseI=TRUE, 
           PhaseII=Phase.Latest %in% c("Phase II Clinical Trial", "Phase III Clinical Trial", "Pre-registration", "Approved"),
           PhaseIII=Phase.Latest %in% c("Phase III Clinical Trial", "Pre-registration", "Approved"),
           Approved=lApprovedUS.EU)
  trial_pairs <- cbind(c("PhaseIN", "PhaseIIN", "PhaseIIIN"),
                       c("PhaseII", "PhaseIII", "Approved"))
  sources <- c("gwas", "omim", "any")
  df <- vector("list", length(sources))
  #risk ratio tests.
  for (i in 1:length(sources)) {
    source_df <- vector("list", length(trial_pairs))
    for (j in 1:nrow(trial_pairs)) {
      # something (hopefully) going wrong at this point.
      source_df[[j]] <- phase_transition_relative_risk_df(annotated_candidates_to_advance, trial_pairs[j,1], trial_pairs[j,2], sources[i], similarity_cutoff, ngene_cutoff)
    }
    df[[i]] <- bind_rows(source_df)
  }
  df <- bind_rows(df)
  df$FirstPhase <- gsub("N", "", gsub("PhaseI", "Phase I", df$FirstPhase))
  df$SecondPhase <- gsub("PhaseI", "Phase I", df$SecondPhase)
  df <- mutate(df, Source=source_names[Source])
  return(bind_rows(df))
}

pairwise_top_status <- function(status_vec1, status_vec2) {
  if (!all(c(status_vec1, status_vec2) %in% c("Approved", "Pre-registration", "Phase I Clinical Trial", "Phase II Clinical Trial", "Phase III Clinical Trial"))) {
    stop("Unsupported status")
  }
  pairwise_top_status_vec <- rep("Phase I Clinical Trial", length(status_vec1))
  pairwise_top_status_vec[status_vec2=="Phase II Clinical Trial"|status_vec1=="Phase II Clinical Trial"] <- "Phase II Clinical Trial" 
  pairwise_top_status_vec[status_vec2=="Phase III Clinical Trial"|status_vec1=="Phase III Clinical Trial"] <- "Phase III Clinical Trial" 
  pairwise_top_status_vec[status_vec2=="Pre-registration"|status_vec1=="Pre-registration"] <- "Pre-registration"
  pairwise_top_status_vec[status_vec2=="Approved"|status_vec1=="Approved"] <- "Approved"
  return(pairwise_top_status_vec)
}

###########################################################
##### Figures #############################################
###########################################################

# Functions to make the main figures of the paper.
# Tend to have limited generalization, but
# Done so we can more easily recreate with, say, another similarity cutoff

make_collected_results_plot <- function(table1_list, 
                                        analysis_classification_vec = c("Updated Clinical"="Mix",
                                                                        "New Clinical"="New",
                                                                        "2013 Data"="Old",
                                                                        "New Genetic"="New",
                                                                        "Clinical Progression"="Pro",
                                                                        "Updated Genetic"="Mix")) {
  for (i in 1:length(table1_list)) {
    table1_list[[i]] <- data.frame(table1_list[[i]], Analysis = names(table1_list)[i], stringsAsFactors = FALSE)
  }
  collected_results <- bind_rows(table1_list)
  
  collected_results <- collected_results %>% mutate(Transition=paste(FirstPhase, SecondPhase, sep=" to\n")) %>% 
    filter(!grepl("&", Source), !FirstPhase=="Preclinical") %>% mutate(Source=ifelse(grepl("GWAS", Source), "GWAS", Source))
  
  collected_results$Transition <- factor(collected_results$Transition, levels=c("Phase I to\nPhase II",   "Phase II to\nPhase III", "Phase III to\nApproved", "Phase I to\nPhase III",  "Phase I to\nApproved"))
  
  collected_results <- collected_results %>% mutate(Type=analysis_classification_vec[Analysis])
  
  #Redone version such that old results are separate.
  collected_new_results <- filter(collected_results, !Type=="Old") %>% dplyr::select(Transition, Source, Estimate, Lower, Upper, Analysis)
  collected_old_results <- filter(collected_results, Type=="Old") %>% dplyr::select(Transition, Source, Estimate, Lower, Upper)
  collected_new_results <- collected_new_results %>% left_join(collected_old_results, by=c("Transition", "Source"), suffix=c("", "Old"))
  
  collected_new_results$Analysis <- factor(collected_new_results$Analysis, 
                                           levels = names(table1_list)[!analysis_classification_vec[names(table1_list)]=="Old"])
  
  combo_plot_scatter2 <- qplot(data=collected_new_results, x=EstimateOld, y=Estimate, ylab="Current Progression Risk Ratio", xlab="Nelson et al. 2015 Progression Risk Ratio", color=Transition, shape=Transition) + 
    facet_grid(Source~Analysis) + 
    geom_hline(yintercept = 1, linetype="dashed") + scale_y_log10(breaks=c(0.25, 0.5,1,2,4)) + 
    scale_x_log10(breaks=c(0.25, 0.5,1,2,4)) + geom_abline(slope=1, intercept=0) +
    geom_errorbarh(aes(xmin=LowerOld, xmax=UpperOld)) +
    geom_errorbar(aes(ymin=Lower, ymax=Upper)) + 
    theme(legend.position = "bottom") + guides(color=guide_legend(ncol=3)) + background_grid(major = "xy", minor = "none")
  
  g <- ggplot_gtable(ggplot_build(combo_plot_scatter2))
  stript <- which(grepl('strip-t', g$layout$name))
  stripr <- which(grepl('strip-r', g$layout$name))
  fillst <- c("#619CFF","#00BA38", "#619CFF", "white")
  fillsr <- c("lightgray", "lightgray")
  k <- 1
  for (i in stript) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fillst[k]
    k <- k+1
  }
  k <- 1
  for (i in stripr) {
    j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
    g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- fillsr[k]
    k <- k+1
  }
  return(g)
}

###########################################################
##### Modeling ############################################
###########################################################

# Portion of design matrix controlling gene level covariates.
gene_design_matrix <- function(target_indication_table, gene_table, p=integer(0), norm=NULL) {
  if (!nrow(gene_table)==length(unique(gene_table$ensembl_id))) {
    stop("Expect one row per gene in gene_table")
  }
  # For those not appearing in p, just join them to the table
  joined_matrix <- left_join(target_indication_table %>% 
                              dplyr::select(ensembl_id), gene_table, by="ensembl_id") %>% 
    dplyr::select(-ensembl_id) %>% as.matrix()
  # For those with names in p, create centered and scaled covariates to order p.
  new_norm <- initialize_norm(norm, length(p))
  for (i in seq_along(p)) {
    cur_col <- which(colnames(joined_matrix)==names(p)[i])
    if (length(cur_col)==1) {
      if (!is.null(norm)) {
        norm_row <- filter(norm, var==names(p)[i]) # will these stably be the same?
        new_vals <- (joined_matrix[,cur_col] - norm_row$mx) / norm_row$sx
        new_cov <- order_p_design_matrix(new_vals, p[i], names(p)[i], standardize.x = FALSE)
        new_norm[[i]] <- norm_row
      } else {
        new_cov <- order_p_design_matrix(joined_matrix[,cur_col], p[i], names(p)[i], standardize.x = TRUE)   
        new_norm[[i]] <- new_cov$norm
      }
      joined_matrix <- cbind(joined_matrix[,-cur_col,drop=FALSE], new_cov$X)
    }
  }
  return(list(X=joined_matrix, norm=bind_rows(new_norm)))
}

# Should support no genetic evidence
# Do this by nothing in filter_list
# If norm is input, it is assumed to be a dataframe.

genetic_evidence_design_matrix <- function(target_indication_table, association_table, MSH_similarity, filter_list, p, norm=NULL) {
  annotated_target_indication_table <- annotate_target_indication_table_with_genetic_evidence(target_indication_table, 
                                                                                              association_table, 
                                                                                              MSH_similarity, 
                                                                                              gene_col_name = "ensembl_id",
                                                                                              filter_list = filter_list)
  
  new_norm <- initialize_norm(norm, length(filter_list))
  if (length(filter_list) > 0) {
    if (length(p)==1) {
      p <- rep(p, length(filter_list))
    }
    assoc_classes <- names(filter_list)
    X_assoc_class <- matrix(nrow = nrow(annotated_target_indication_table), ncol = length(assoc_classes))
    colnames(X_assoc_class) <- gene_intercept_colname(assoc_classes)
    X_similarity <- matrix(nrow = nrow(annotated_target_indication_table), ncol = 0)
    for (i in 1:length(assoc_classes)) {
      mesh <- annotated_target_indication_table[[similarity_base(assoc_classes[i])]]
      gene <- annotated_target_indication_table[[gene_intercept_colname(assoc_classes[i])]]
      # mesh standardization special: standardize by mean and sd of nonzero values
      if (p[i] > 0) {
        if (is.null(norm)) {
          mesh <- standardize_mesh(mesh, gene)
          new_norm[[i]] <- data.frame(var=similarity_base(assoc_classes[i]), mx=mesh$mx, sx=mesh$sx, stringsAsFactors = FALSE)
        } else {
          norm_row <- filter(norm, var==similarity_base(assoc_classes[i]))
          mesh <- list(res = gene * (mesh - norm_row$mx) / norm_row$sx, mx = norm_row$mx, sx = norm_row$sx)
          new_norm[[i]] <- norm_row
        }
        mesh_mat <- order_p_design_matrix(x = mesh$res, p = p[i], name_base = similarity_base(assoc_classes[i]), standardize.x = FALSE)
        X_similarity <- cbind(X_similarity, mesh_mat$X)
      }
      X_assoc_class[,gene_intercept_colname(assoc_classes[i])] <- as.integer(gene)
    }
    X <- cbind(X_assoc_class, X_similarity)
  } else {
    X <- matrix(nrow=nrow(annotated_target_indication_table), ncol=0)
  }
  
  return(list(X=X, norm=bind_rows(new_norm)))
}

initialize_norm <- function(norm, n) {
  new_norm <- vector("list", n)
  return(new_norm)
}

order_p_design_matrix <- function(x, p, name_base, standardize.x=TRUE) {
  n <- length(x)
  if (standardize.x) {
    x <- center_scale_x(x)
  } else {
    x <- list(res=x, mx=0, sx=1)
  }
  X <- matrix(nrow=n, ncol=p)
  for (i in seq_len(p)) {
    X[,i] <- x$res^i
  }
  if (p > 0) {
    colnames(X) <- paste0(name_base, 1:p)
  }
  return(list(X=X, norm=data.frame(var = name_base, mx = x$mx, sx = x$sx, stringsAsFactors = FALSE)))
}

center_scale_x <- function(x) {
  mx <- mean(x, na.rm=TRUE)
  sx <- sd(x, na.rm=TRUE)
  z <- (x-mx)/sx
  return(list(res = z, mx = mx, sx = sx))
}

standardize_mesh <- function(mesh, gene, vals_to_standardize=mesh, genes_to_standardize=gene) {
  mean_mesh <- mean(mesh[gene])  
  sd_mesh <- sd(mesh[gene])
  std <- (vals_to_standardize - mean_mesh) / sd_mesh
  std[!genes_to_standardize] <- 0
  return(list(res= std, mx=mean_mesh, sx=sd_mesh))
}

# Function to construct a design matrix indicating whether or not 
# each human disease top MeSH term applies to
# each mesh term in target-indication table

top_mesh_design_matrix <- function(target_indication_table, top_mesh) {
  # Exclude nondisease, animal disease (C22)
  utop <- unique(top_mesh$MSH.Top)
  X_msh <- matrix(nrow = nrow(target_indication_table), ncol=length(utop), data =0)
  colnames(X_msh) <- utop
  for (i in seq_along(utop)) {
    X_msh[target_indication_table$MSH %in% (filter(top_mesh, MSH.Top==utop[i]) %>% pull(Name)), i] <- 1
  }
  colnames(X_msh) <- make.names(colnames(X_msh))
  return(list(X=X_msh, norm=tibble()))
}

# Function to run stan from input data frames.  Returns a list with components stan_out and stan_in.
# I am electing not to store data now as I am not sure what can be shared, but we need normalization 
# in order to predict new observations.
# ... is additional arguments to sampler.  Optionally, save output to file.

run_stan <- function(target_indication_table, association_table, MSH_similarity,
                     top_mesh, gene_table, pg, pe, filter_list, 
                     sigmab, mualpha, sigmaalpha, save_file=NULL, stan_path="../src/LogisticRegression.stan",...) {
  data = list(target_indication_table=target_indication_table, association_table=association_table, 
              MSH_similarity=MSH_similarity, top_mesh=top_mesh, gene_table=gene_table)
  params = list(pg=pg, pe=pe, filter_list=filter_list, sigmab=sigmab, sigmaalpha=sigmaalpha, mualpha=mualpha)
  stan_input <- do.call(generate_stan_input, c(data, params))
  stan_out <- stan(file=stan_path, data=stan_input$input,...)
  res <- list(stan_input=stan_input$input, stan_out=stan_out, params=params, norm=stan_input$norm)
  if (!is.null(save_file)) {
    saveRDS(res, file = save_file)
  }
  return(res)
}

# Taking the inputs to functions for generating mesh, gene, and genetic evidence design matrices
# and prior parameters and combining resulting outputs to get design matrix.

generate_stan_input <- function(target_indication_table, association_table, MSH_similarity,
                                top_mesh, gene_table, pg, pe, filter_list, 
                                sigmab, mualpha, sigmaalpha, norm=NULL) {
  dmat <- generate_full_design_matrix(target_indication_table, association_table, MSH_similarity,
                                      top_mesh, gene_table, pg, pe, filter_list, norm)
  missing_indices <- apply(is.na(dmat$X), 1, any)
  y <- as.integer(target_indication_table$lApprovedUS.EU)
  X_nonmissing <- dmat$X[!missing_indices,,drop=FALSE]
  y_nonmissing <- y[!missing_indices]
  stan_input <- list(X=X_nonmissing, N=nrow(X_nonmissing), d=ncol(X_nonmissing), y=y_nonmissing, 
                     sigmab=sigmab, mualpha=mualpha, sigmaalpha=sigmaalpha)
  return(list(input = stan_input, norm=dmat$norm))
}

generate_full_design_matrix <- function(target_indication_table, association_table, MSH_similarity,
                                        top_mesh, gene_table, pg, pe, filter_list, norm=NULL) {
  X_msh <- top_mesh_design_matrix(target_indication_table, top_mesh)
  X_gene <- gene_design_matrix(target_indication_table, gene_table, p=pg, norm=norm)
  X_ge <- genetic_evidence_design_matrix(target_indication_table, association_table, MSH_similarity, filter_list, pe, norm=norm)
  X <- cbind(X_ge$X, X_msh$X, X_gene$X)
  norm <- bind_rows(list(X_msh$norm, X_gene$norm, X_ge$norm))
  return(list(X=X, norm=norm))
}

# get the odds ratio for drugs with genetic evidence versus with no genetic association for gene
# for different values of trait-gene similarity

generate_prediction_df <- function(stan_input, stan_out, params, norm, similarities) {
  if (length(params$pe)==1) {
    params$pe <- rep(params$pe, length(params$filter_list))
  }
  pred_list <- vector("list", length(params$filter_list))
  for (i in 1:length(params$filter_list)) {
    gene_col <- which(colnames(stan_input$X) == gene_intercept_colname(names(params$filter_list)[i]))
    mesh_base <- similarity_base(names(params$filter_list)[i])
    norm_row <- filter(norm, var==mesh_base)
    similarity_std <- (similarities - norm_row$mx) / (norm_row$sx)
    posterior_samps_intercept <- as.matrix(stan_out, pars=paste0("beta[", gene_col, "]"))
    predicted_effect <- matrix(data = posterior_samps_intercept[,1], nrow=nrow(posterior_samps_intercept), 
                               ncol=length(similarities), byrow=FALSE)
    for (j in seq_len(params$pe[i])) {
      mesh_col <- which(colnames(stan_input$X)==paste0(similarity_base(names(params$filter_list)[i]), j))
      posterior_samps <- as.matrix(stan_out, 
                                   pars=paste0("beta[", mesh_col, "]"))
      predicted_effect <- predicted_effect + posterior_samps[,1,drop=FALSE] %*% matrix(data=similarity_std^j, nrow=1)
    }
    pred_list[[i]] <- data.frame(Median = exp(apply(predicted_effect, 2, median)), Mean = apply(exp(predicted_effect), 2, mean),
                                 Lower = exp(apply(predicted_effect, 2, quantile, 0.025)), Upper = exp(apply(predicted_effect, 2, quantile, 0.975)),
                                 Source = names(params$filter_list)[i], Similarity = similarities, stringsAsFactors = FALSE)
  }
  return(bind_rows(pred_list))
}

# Function that takes an ensembl id and a MeSH term and returns a predicted success probability.  Requires input and output from a stan run, and 
# confidence level, as well as some data tables allowing to get covariates for the genes.
# TODO: new_dev_time not currently used.
# TODO: we could probably speed it up alot using changes to similarity function such that it only computes the parts actually needed (1 row here)


# parameter cutoff determines the minimum similarity for which ensembl ids and mesh term pairs are included.
predict_gene_mesh <- function(ensembl_ids, mesh_terms, conf, stan_res, similarity, gene_table, top_mesh, association_table, default_time = 13903.5, cutoff=0.5, find_similar_trait=TRUE) {
  X_table <- new_design_rows_df(ensembl_id = ensembl_ids, mesh_term = mesh_terms, association_table = association_table, 
                                 MSH_similarity = similarity, top_mesh = top_mesh, gene_table = gene_table, stan_res = stan_res, cutoff=cutoff)
  X <- as.matrix(X_table[,!colnames(X_table) %in% c("ensembl_id", "MSH")])
  posterior_samps_beta <- as.matrix(stan_res$stan_out, pars=paste0("beta[", 1:ncol(X), "]"))
  posterior_samps_alpha <- matrix(data = as.matrix(stan_res$stan_out, pars="alpha"), nrow = nrow(posterior_samps_beta), ncol = nrow(X), byrow=FALSE)
  pred <- posterior_samps_alpha + posterior_samps_beta %*% t(X)
  pred_p <- 1 / (1 + exp(-pred))
  pres <- data.frame(mean = apply(pred_p, 2, mean), median = apply(pred_p, 2, median), 
                     lower = apply(pred_p, 2, quantile, (1-conf)/2), upper = apply(pred_p, 2, quantile, 1-(1-conf)/2))
  # I would also like the component of beta specifically due to genetic evidence.g
  if (find_similar_trait) {
    annotated_table <- annotate_target_indication_table_with_genetic_evidence(X_table[,c("ensembl_id", "MSH")], 
                                                                              association_table = association_table,
                                                                              MSH_similarity = similarity, 
                                                                              gene_col_name = "ensembl_id")
    trait_res <- vector("list", length(stan_res$params$filter_list))
  }
  sim_res <- vector("list", length(stan_res$params$filter_list))
  ge_col <- c()
  for (i in seq_along(stan_res$params$filter_list)) {
    assoc_class <- names(stan_res$params$filter_list)[i]
    ge_col <- c(ge_col, which(colnames(X)==gene_intercept_colname(assoc_class)))
    # TODO: absorb into an input checking function OR absorb into function to get column names?
    if (length(stan_res$params$pe)==1) {
      stan_res$params$pe <- rep(stan_res$params$pe, length(stan_res$params$filter_list))
    }
    if (stan_res$params$pe[i] > 0) {
      sim_col <- similarity_base(assoc_class)
      ge_col <- c(ge_col, which(colnames(X) %in% paste0(sim_col, 1:stan_res$params$pe[i])))
      sims <- ifelse(X[,gene_intercept_colname(assoc_class)]==1, X[,paste0(sim_col, 1)] * filter(stan_res$norm, var==sim_col)$sx + filter(stan_res$norm, var==sim_col)$mx, NA)
    } else {
      sims <- NA
    }
    sim_res[[i]] <- data.frame(Similarity = sims, Source = rep(toupper(assoc_class), length(sims)), ensembl_id = X_table$ensembl_id,
                               MSH = X_table$MSH, stringsAsFactors = FALSE)
    if (find_similar_trait) {
      trait <- annotated_table[[trait_colname(assoc_class)]]
      trait_res[[i]] <- data.frame(ensembl_id = X_table$ensembl_id, MSH = X_table$MSH, Source = rep(toupper(assoc_class), nrow(annotated_table)),
                                   Trait = trait, stringsAsFactors = FALSE)
    }
  }
  ge_component <- exp(posterior_samps_beta[,ge_col] %*% t(X[,ge_col]))
  sim_res <- bind_rows(sim_res) %>% spread(Source, Similarity)
  stan_res <- data.frame(p = pres, OR = apply(ge_component, 2, mean), OR.lower = apply(ge_component, 2, quantile, (1-conf)/2), OR.upper = apply(ge_component, 2, quantile, 1-(1-conf)/2), X_table[,c('ensembl_id', 'MSH')], stringsAsFactors = FALSE)
  res <- left_join(stan_res, sim_res)
  if (find_similar_trait) {
    trait_res <- bind_rows(trait_res) %>% spread(Source, Trait)
    res <- left_join(res, trait_res, by = c("ensembl_id", "MSH"), suffix = c("Similarity", "Trait"))
  }
  return(res)
}

# Time is optionally set to a default value (default is do this).  This is recommended for new targets in particular as a certain amount of
# time is required for a target to have a positive probability of success.  Set to NA if you don't want.
# Default time is 10 before dataset download. 
# What this means as we estimate the success for a target for which Pharmaprojects knows about 10 years of development.
# This is a commonly stated figure and consistent with calculations from Pharmaprojects 
# (calculating time between drug with target first added and first launched or registered, excluding negative values, get 10.35 years on average).

new_design_rows_df <- function(ensembl_ids, mesh_terms, association_table, MSH_similarity, gene_table, top_mesh, stan_res, 
                           default_time = 13903.5, cutoff=0.5) {
  tmp_table <- select_gene_msh_pairs(ensembl_ids, mesh_terms, association_table, MSH_similarity, cutoff) %>% mutate(lApprovedUS.EU=NA)
  new_rows <- matrix(data=0, ncol=ncol(stan_res$stan_input$X), nrow=nrow(tmp_table))
  colnames(new_rows) <- colnames(stan_res$stan_input$X)
  # generate the input again using original normalization.
  X <- generate_full_design_matrix(tmp_table, association_table = association_table, MSH_similarity = MSH_similarity, 
                                   gene_table = gene_table, top_mesh = top_mesh,
                                   pg = stan_res$params$pg, pe = stan_res$params$pe, filter_list = stan_res$params$filter_list, norm=stan_res$norm)$X
  # Any ones not in there are (currently) MeSH headings not found in the row and can be zero
  # Any ones not in the model are excluded
  new_rows[,intersect(colnames(X), colnames(new_rows))] <- X[,intersect(colnames(X), colnames(new_rows))]
  # Any ones that are NA are currently set to zero, equivalent to sample mean for continuous variables.
  # If there is a default, set to the standardized version of that.
  na_ind <- which(is.na(new_rows), arr.ind=TRUE)
  new_rows[na_ind] <- 0
  if (!is.na(default_time)) {
    default_time <- stan_res$norm %>% filter(var=="Time") %>% mutate(default = (default_time - mx) / sx) %>% pull(default)
    if (length(stan_res$params$pg["Time"])==0) {
      new_rows[,"Time"] <- default_time     
    } else {
      for (i in 1:stan_res$params$pg["Time"]) {
        new_rows[,paste0("Time", i)] <- default_time^i
      } 
    }
  }
  new_rows_df <- data.frame(tmp_table[,1:2], new_rows, stringsAsFactors = FALSE)
  return(new_rows_df)
}

select_gene_msh_pairs <- function(ensembl_ids, mesh_terms, association_table, MSH_similarity, cutoff) {
  cand_ids <- ensembl_ids[ensembl_ids %in% association_table$ensembl_id]
  cand_pairs <- vector("list", length(cand_ids))
  for (i in seq_along(cand_ids)) {
    associated_MSH <- filter(association_table, ensembl_id == cand_ids[i]) %>% pull(MSH) %>% unique()
    msh_cols <- MSH_similarity[mesh_terms,associated_MSH,drop=FALSE]
    similar_mesh <- rownames(msh_cols)[apply(msh_cols >= cutoff, 1, any)]
    if (length(similar_mesh) > 0) {
      cand_pairs[[i]] <- data.frame(ensembl_id = cand_ids[i], MSH = similar_mesh, stringsAsFactors = FALSE)
    }
  }
  if (all(sapply(cand_pairs, is.null))) {
    cand_pairs <- data.frame(ensembl_id = character(0), MSH = character(0), stringsAsFactors = FALSE)
  } else {
    cand_pairs <- bind_rows(cand_pairs)
  }
  return(cand_pairs)
}
