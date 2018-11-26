#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(DT)

path = "../"

stan_res <- readRDS(paste0(path, "results/ShinyAppPrecomputed.rds"))
gene_trait_association_table <- 
  read.delim(paste0(path, "data/gene_trait_assoc.tsv"), stringsAsFactors = FALSE)
TI_table <- 
  read.delim(paste0(path, "data/target_indication.tsv"), stringsAsFactors = FALSE)
stan_res <- stan_res %>% 
  left_join(gene_trait_association_table %>% 
              dplyr::select(ensembl_id, symbol) %>% distinct) %>% 
  filter(!is.na(symbol)) %>% 
  mutate(PharmaprojectsTarget=ensembl_id %in% TI_table$ensembl_id)

source(paste0(path, "/src/MeSHFunctions.R"))
source(paste0(path, "/src/StatisticalFunctions.R"))

ui <- fluidPage(
   
   # Application title
   titlePanel("Drug Target Success Predictions"),
   
   tabsetPanel(
     tabPanel("Genetic Evidence by Gene", fluid = TRUE,
              sidebarLayout(
                sidebarPanel(
#                  selectInput("vocab", "Gene name format", tolower(keytypes(org.Hs.eg.db)), selected = "symbol", multiple = FALSE,
#                              selectize = TRUE, width = NULL, size = NULL), 
                  selectizeInput("target",
                                "Target",
                                unique(stan_res$symbol), 
                                selected = "IL10",
                                multiple = FALSE, width = NULL, size = NULL),
                  selectizeInput("model", 
                                 "Model", 
                                 c("GWAS and OMIM", "GWAS only", "OMIM only"), 
                                 selected = "GWAS and OMIM", 
                                 multiple = FALSE, width = NULL, size = NULL),
                  width=2
                ),
                mainPanel(
                  DTOutput("PredictionSummaryGene")
                )
              )
     ),
     tabPanel("Genetic Evidence by Indication", fluid = TRUE,
              sidebarLayout(
                sidebarPanel(
                  selectizeInput("msh", "Indication (2017 MeSH term)", unique(stan_res$MSH), selected = "Psoriasis", multiple = FALSE,
                              width = NULL, size = NULL), 
                  selectizeInput("model2", "Model", c("GWAS and OMIM", "GWAS only", "OMIM only"), selected = "GWAS and OMIM", multiple = FALSE,
                              width = NULL, size = NULL),
#                  selectInput("vocab2", "Gene name format", tolower(keytypes(org.Hs.eg.db)), selected = "symbol", multiple = FALSE,
#                              selectize = TRUE, width = NULL, size = NULL), 
                  checkboxInput("knowntargets", "Exclude Non-targets", value = FALSE, width = NULL),
                  width=2
                ),
                mainPanel(
                  DTOutput("PredictionSummaryIndication")
                )
              )
     )
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  output$PredictionSummaryGene <- renderDT({
    selected_model <- c("GWAS and OMIM"="both", "GWAS only"="gwas", "OMIM only"="omim")[input$model]
    res <- filter(stan_res, symbol==input$target, Model==selected_model) %>% 
      dplyr::select(MSH, p.mean, OR, GWASSimilarity, GWASTrait, OMIMSimilarity, OMIMTrait) %>% 
      arrange(desc(OR))
    datatable(res, colnames = c('GWAS Trait' = 'GWASTrait', 'OMIM Trait' = 'OMIMTrait', 
                                     'GWAS Similarity' = 'GWASSimilarity', 'OMIM Similarity' = 'OMIMSimilarity', 
                                     'Success Probability' = 'p.mean', 'Genetic Evidence (Odds Ratio)' = 'OR', 
                                'Indication (MeSH)' = 'MSH')) %>%
      formatRound('Success Probability',2) %>%
      formatRound('Genetic Evidence (Odds Ratio)', 2) %>%
      formatRound('GWAS Similarity', 2) %>%
      formatRound('OMIM Similarity', 2)
   })
  output$PredictionSummaryIndication <- renderDT({
    selected_model2 <- c("GWAS and OMIM"="both", "GWAS only"="gwas", "OMIM only"="omim")[input$model2]
    res <- filter(stan_res, MSH==input$msh, Model==selected_model2) %>% 
      dplyr::select(symbol, p.mean, OR, GWASSimilarity, GWASTrait, OMIMSimilarity, OMIMTrait, PharmaprojectsTarget) %>% 
      arrange(desc(OR))
    if (input$knowntargets) {
      res <- filter(res, PharmaprojectsTarget)
    }
    datatable(res, colnames = c('GWAS Trait' = 'GWASTrait', 'OMIM Trait' = 'OMIMTrait', 
                                'GWAS Similarity' = 'GWASSimilarity', 'OMIM Similarity' = 'OMIMSimilarity', 
                                'Success Probability' = 'p.mean', 'Genetic Evidence (Odds Ratio)' = 'OR', 
                                'Target' = 'symbol', 'Pharmaprojects Target' = 'PharmaprojectsTarget')) %>%
      formatRound('Success Probability',2) %>%
      formatRound('Genetic Evidence (Odds Ratio)', 2) %>%
      formatRound('GWAS Similarity', 2) %>%
      formatRound('OMIM Similarity', 2)
  })
}

# Run the application 
shinyApp(ui = ui, server = server)

