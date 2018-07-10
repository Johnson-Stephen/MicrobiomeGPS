# 
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
#
#
# Adding more comments -- test
library(MASS) 
library(ggbiplot)
library(ape)
library(vegan)
library(pscl)
library(glmmADMB)
library(aod)
library(nlme)
library(MiRKAT)
library(matrixStats)
library(gplots)
library(scales)
library(ggplot2)
library(GUniFrac)
library(rpart)
library(qvalue)
library(phangorn)
library(RColorBrewer)
library(squash)
library(rhdf5)
library(biom)
library(randomForest)
library(Boruta)
library(ade4)
library(Daim)
library(phyloseq)
library(xtable)
library(shiny)
library(shinydashboard)
library(Matrix)
library(ggnet)
library(igraph)
library(gridExtra)
library(devtools)
library(SpiecEasi)
library(ROCR)
library(cluster)
library(dplyr)
library(data.table)
library(tidyverse)
library(biom)
library(phyloseq)
library(DT)
library(iheatmapr)
library(RColorBrewer)
library(reshape2)
library(plotly)
options(shiny.maxRequestSize = 100*1024^2)
shinyApp(
  ui = dashboardPage(
    dashboardHeader(title = "MicrobiomeGPS"),
    dashboardSidebar(
      sidebarMenu(
        menuItem("Create Dataset", tabName = "create_dataset", icon = icon("dashboard")),
        menuItem("Summary statistics", tabName = "summary_statistics", icon = icon("th")),
        menuItem("Alpha diversity", tabName = "alpha_diversity", icon = icon("th")),
        menuItem("Beta diversity", tabName = "beta_diversity", icon = icon("th")),
        menuItem("Taxa diversity", tabName = "taxa_diversity", icon = icon("th")),
        menuItem("Predictive modeling", tabName = "predictive_modeling", icon = icon("th")),
        menuItem("Functional analysis", tabName = "functional_analysis", icon = icon("th")),
        menuItem("Subtype discovery", tabName = "subtype_discovery", icon = icon("th")),
        menuItem("Network analysis", tabName = "network_analysis", icon = icon("th")),
        menuItem("Generate reports", tabName = "Report", icon = icon("th"))
      )
    ),
    dashboardBody(#"MicrobiomeGPS",
      height="100%",
      tabItems(
        tabItem(tabName = "create_dataset",
                fluidRow(
                  column(width=4,
                         box(title="1. Upload files", solidHeader=TRUE, status="primary", width=NULL,
                             #h2("1. Upload files"),
                             fileInput("mapping_file", label = h3("Mapping file")),
                             fileInput("tree_file", label = h3("Tree file")),
                             fileInput("biom_file", label = h3("Biom file")),
                             fileInput("kegg_file", label = h3("KEGG file")),
                             fileInput("cog_file", label = h3("COG file"))
                         ),
                         box(title="2. Filter data", solidHeader=TRUE, status="primary", width=NULL,
                             numericInput("filter_dep", "Remove samples with sequence count below:", 2000, min = 0, step = 1000),
                             numericInput("otu_filter_count", "Remove OTUs with sequence count at or below:", 1, min = 0, step = 100),
                             numericInput("otu_filter_freq", "Remove OTUs with occurence frequency at or below:", 0, min = 0, max = 100, step = 1)
                         ),
                         box(title="3. Normalize data", solidHeader=TRUE, status="primary", width=NULL,
                             h4("(I) Rarefied object:"),
                             numericInput("rare_dep", "Rarefication depth:", 10000, min = 0, step = 100),
                             h4("(II) Unrarefied object"),
                             selectInput("size_factor", "Size Factor:", choices=c("TSS", "CSS", "GMPR", "RLE", "TMM"), selected=c("TSS"))
                         )
                  ),
                  column(width=4,
                         box(title="4. Subset data", solidHeader=TRUE, status="primary", width=NULL,
                             uiOutput("maptest"),
                             uiOutput("maptest2"),
                             uiOutput("maptest6"),
                             uiOutput("maptest3"),
                             tags$style(type="text/css", "textarea {width:100%}"),
                             tags$label('for'='input_text', "List of samples to exclude (Case sensitive):"),
                             tags$textarea(id = 'input_text', placeholder = 'Paste Sample IDs here', rows = 6, ""),
                             tags$label('for'='output_text', "Invalid sample names, check input:"),
                             verbatimTextOutput("output_text"),
                             tags$style(type="text/css", "textarea {width:100%}"),
                             tags$label('for'='r_select', "Enter R expression to select samples (e.g., 'Sex == 'M' & BMI < 50')"),
                             tags$textarea(id = 'r_select', placeholder = 'Enter R expression here', rows = 6, NULL),
                             tags$label('for'='r_select_output', "Invalid sample names, check input:"),
                             verbatimTextOutput("r_select_output")
                         ),
                         box(title="5. Create dataset", solidHeader=TRUE, status="primary", width=NULL,
                             actionButton("create_dataset", label = "Submit"),
                             h3(textOutput("status", container = span))
                             
                         )
                         
                  ),
                  column(width=4,
                         box(title="6. Assign variable types", solidHeader=TRUE, status="primary", width=NULL,
                             uiOutput("select_cat_vars"),
                             uiOutput("select_con_vars"),
                             uiOutput("select_ord_vars")
                         ),
                         box(title="7. Select variables", solidHeader=TRUE, status="primary", width=NULL,
                             uiOutput("select_var"),
                             uiOutput("select_covars"),
                             uiOutput("select_subject")
                         )
                  )
                )
        ),
        tabItem(tabName = "summary_statistics",
                fluidRow(
                  box(width=3,
                      numericInput("prev", "Minimum prevalence threshold (%):", 10, min = 0, max = 100, step = 10),
                      numericInput("abund", "Minimum abundance threshold (%):", 0.2, min = 0, max = 100, step = 1),
                      actionButton("summary_stats", label = "Submit"),
                      h3(textOutput("summary_status", container = span))
                  ),
                  box(width=9,
                      tabBox(width=9,
                        tabPanel("Sequencing statistics",
                                 plotOutput("cov_dist", width=900, height=600),
                                 plotOutput("cov_boxplot", width=900, height=600)
                        ),
                        tabPanel("Prevalence",
                                 DT::dataTableOutput("phy.prev", width="100%"),
                                 DT::dataTableOutput("fam.prev", width="100%"),
                                 DT::dataTableOutput("gen.prev", width="100%")
                        ),
                        tabPanel("Abundance",
                                 DT::dataTableOutput("phy.abund", width="100%"),
                                 DT::dataTableOutput("fam.abund", width="100%"),
                                 DT::dataTableOutput("gen.abund", width="100%")
                        ),
                        tabPanel("Barplots",
                                 selectInput("barplot_level", "Level:", choice=c("Phylum", "Family", "Genus"), selected='Phylum'),
                                 plotOutput("summary_barplot1"),
                                 plotOutput("summary_barplot2")
                        ),
                        tabPanel("Heatmaps",
                                 selectInput("heatmap_type", "Heatmap Type", choice=c("Proportional", "Binary", "Ranked"), selected='Proportional'),
                                 iheatmaprOutput("summary_heatmap")
                        )
                      )
                  )
                )
        ),
        tabItem(tabName="alpha_diversity",
                fluidRow(
                  box(width=3,
                      h2("Alpha diversity"),
                      uiOutput("alpha_measures"),
                      numericInput("alpha_rare_dep", "Rarefication depth:", 10000, min = 0, step = 100),
                      numericInput("rare_iter", "Rarefication iterations:", 5, min = 1, step = 1),
                      selectInput("alpha_nonrare_test", label = "Nonrarefied Test", choices = c("NO","YES")),
                      actionButton("run_alpha", label="Submit"),
                      verbatimTextOutput("alpha_text")
                  ),
                  box(width=9,
                      tabBox(width=9,
                             tabPanel("Rarefaction",
                                      fluidRow(
                                        h2("Rarefaction Curve"),
                                        plotOutput("rarefy_curve", width=900, height=600),
                                        h2("Rarefaction Boxplot"),
                                        plotOutput("rarefy_boxplot", width=900, height=600)
                                      )
                             ),
                             tabPanel("Association test",
                                      uiOutput("alpha_association_tab1"),
                                      uiOutput("alpha_association_tab2"),
                                      uiOutput("alpha_association_tab3"),
                                      uiOutput("alpha_association_tab4"),
                                      uiOutput("alpha_association_tab5"),
                                      uiOutput("alpha_association_tab6"),
                                      uiOutput("alpha_association_tab7"),
                                      uiOutput("alpha_association_tab8")
                             )
                      )
                  )
                )
        ),
        tabItem(tabName="beta_diversity",
                fluidRow(
                  box(width=3,
                      h2("Beta diversity"),
                      uiOutput("beta_measures"),
                      checkboxInput("rf_check", "Rarefaction", value = TRUE),
                      selectInput("ord_measure", "Ordination method", choices=c("PCoA", "NDMS"), selected="PCoA"),
                      actionButton("run_beta", label="Submit"),
                      verbatimTextOutput("beta_text")
                  ),
                  box(width=9,
                      tabBox(width=9,
                             tabPanel("Ordination",
                                      plotlyOutput("ordination", width=900, height=600)         
                             ),
                             tabPanel("Distance comparison",
                                      fluidRow(
                                        h2("Distance comparison boxplots"),
                                        plotOutput("distance_comparison_boxplot", width=900, height=600)
                                      )
                             ),
                             tabPanel("Beta Association",
                                      h2("Permanova test"),
                                      p("Permanova is a multivariate analysis of variance based on distance matrices and permutation, partitioning distance matrices among sources of variation and fitting linear models to distance matrices. Permanova analysis was performed using the adonis package in R."),
                                      htmlOutput("permanova"),
                                      h2("MiRKAT test"),
                                      p("MiRKAT is a kernel-based association test based on ecological distance matrices. MiRKAT produces analytic p-values for individual distance metrics, as well as a permutation-based omnibus p-value that combines multiple distance metrics, for a more robust and powerful assessment of significance."),
                                      htmlOutput("mirkat"),
                                      h2("Beta dispersion test"),
                                      p("BETADISPER is part of the R vegan package. It is a multivariate analogue of Levene’s test for homogeneity of varainces. Non-euclidean distances between objects and group centroids are handled by reducing the original distances to principal coordinates."),
                                      htmlOutput("betadisper")
                             )
                      )
                  )
                )
        ),
        tabItem(tabName="taxa_diversity",
                fluidRow(
                  box(width=3,
                      h2("Taxa diversity"),
                      selectInput("taxa_method", "Differential abundance test", choices=c("Permutation" = "perm", "Wilcox" = "wilcox", "Wilcox.pair", "Kruskal-Wallis", "Twopart", "Fisher", "OverdispersedPoisson","OverdispersedBinomial", "NegativeBinomial", "ZeroInflatedNB"), selected="Permutation"),
                      selectInput("normalization", "Normalization method", choices=c("Rarefaction", "TSS", "CSS", "GMPR", "RLE", "TMM"), selected="GMPR"),
                      selectInput("taxa_trans", "Transformation method", choices=c("Square root", "Square root arcsine"), selected="Square root"),
                      selectInput("taxa_outliers", "Addressing outliers", choices=c("Winsorization", "Reweighting"), selected="Winsorization"),
                      numericInput("winsor", "Winsorization quantile", 0.97, min = 0, max = 1.0, step = 0.01),
                      selectInput("mult_test", "Method for multiple testing correction:", choices=c("fdr", "raw", "None", "Bonferroni", "Storey-q"), selected="fdr"),
                      numericInput("sig_level", "Significance level (%)", 10, min = 0, step = 10),
                      numericInput("taxa_prev", "Minimum prevalence threshold (%):", 10, min = 0, max = 100, step = 10),
                      numericInput("taxa_abund", "Minimum abundance threshold (%):", 0.2, min = 0, max = 100, step = 1),
                      uiOutput("vis_level"),
                      actionButton("run_taxa", label="Submit"),
                      verbatimTextOutput("taxa_text")
                  ),
                  box(width=9,
                      tabBox(width=9,
                             tabPanel("Boxplots",
                                      plotOutput("taxa_boxplots")
                             ),
                             tabPanel("Barplots",
                                      plotOutput("taxa_barplots")
                             ),
                             tabPanel("Effect size",
                                      plotOutput("effect_size")
                             ),
                             tabPanel("PCA biplot",
                                      plotOutput("taxa_biplot")
                             ),
                             tabPanel("Heatmaps",
                                      iheatmaprOutput("taxa_prop_heatmap"),
                                      iheatmaprOutput("taxa_rank_heatmap")
                             ),
                             tabPanel("Cladogram",
                                      htmlOutput("cladogram")
                             ),
                             tabPanel("Test results",
                                      htmlOutput("taxa_test_results")
                             )
                      )
                  )
                )
        ),
        tabItem(tabName="predictive_modeling",
                fluidRow(
                  box(width=3,
                      h2("Predictive modeling"),
                      selectInput("pred_method", "Predictive method", choices=c("Random forest"), selected="Random forest"),
                      numericInput("bootstrap_num", "Bootstrap number:", 100, min=0, step = 100),
                      selectInput("boruta_level", "Boruta significance level:", choices=c("Tentative", "Confirmed"), selected="Tentative"),
                      selectInput("pred_level", "Prediction level:", choices=c("OTU", "Genus"), selected="OTU"),
                      selectInput("pred_norm", "Normalization method", choices=c("Rarefaction", "TSS", "CSS", "GMPR", "RLE", "TMM"), selected="GMPR"),
                      selectInput("pred_trans", "Transformation method", choices=c("Square root", "Square root arcsine"), selected="Square root"),
                      selectInput("pred_outliers", "Addressing outliers", choices=c("Winsorization", "Reweighting"), selected="Winsorization"),
                      numericInput("pred_prev", "Minimum prevalence threshold (%):", 10, min = 0, max = 100, step = 10),
                      numericInput("pred_abund", "Minimum abundance threshold (%):", 0.2, min = 0, max = 100, step = 1),
                      actionButton("run_pred", label="Submit"),
                      verbatimTextOutput("pred_text")
                  ),
                  box(width=9,
                      tabBox(width=9,
                             tabPanel("Prediction evaluation",
                                      imageOutput("classification_error", width=900, height=600)
                             ),
                             tabPanel("Bootstrap-validated ROC curves",
                                      imageOutput("bootstrap_roc_genus", width = 900, height = 600),
                                      imageOutput("bootstrap_roc_species", width = 900, height = 600)
                             ),
                             tabPanel("Boruta selected features",
                                      htmlOutput("boruta_features")
                             ),
                             tabPanel("Boruta barplots",
                                      htmlOutput("boruta_barplots_agg"),
                                      htmlOutput("boruta_barplots_ind")
                             ),
                             tabPanel("Boruta boxplots",
                                      htmlOutput("boruta_boxplots")
                             )
                      )
                  )
                )
        ),
        tabItem(tabName="functional_analysis",
                fluidRow(
                  box(width=3,
                      h2("Functional Analysis"),
                      selectInput("func_method", "Functional diversity method", choices=c("perm", "perm.pair", "wilcox", "wilcox.pair", "kruskal", "twopart", "Spearman", "OP", "NB", "ZINB"), selected="perm"),
                      selectInput("func_mult_test", "Method for multiple testing correction:", choices=c("fdr", "raw", "None", "Bonferroni", "Storey-q"), selected="fdr"),
                      numericInput("func_sig_level", "Significance level (%)", 10, min = 0, step = 10),
                      selectInput("func_normalization", "Normalization method", choices=c("Rarefaction", "TSS", "CSS", "GMPR", "RLE", "TMM"), selected="GMPR"),
                      selectInput("func_trans", "Transformation method", choices=c("Square root", "Square root arcsine"), selected="Square root"),
                      selectInput("func_outliers", "Addressing outliers", choices=c("Winsorization", "Reweighting"), selected="Winsorization"),
                      numericInput("func_winsor", "Winsorization quantile", 0.97, min = 0, max = 1.0, step = 0.01),
                      numericInput("func_prev", "Minimum prevalence threshold (%):", 10, min = 0, max = 100, step = 10),
                      numericInput("func_abund", "Minimum abundance threshold (%):", 0.2, min = 0, max = 100, step = 1),
                      #uiOutput("func_vis_level"),
                      actionButton("run_func", label="Submit"),
                      verbatimTextOutput("func_text")
                  ),
                  box(width=9,
                      tabBox(width=9,
                             tabPanel("KEGG barplots",
                                      htmlOutput("kegg_barplot_agg"),
                                      htmlOutput("kegg_barplot_ind")
                             ),
                             tabPanel("KEGG boxplots",
                                      htmlOutput("kegg_boxplot_agg"),
                                      htmlOutput("kegg_boxplot_ind")
                             ),
                             tabPanel("KEGG effect size",
                                      htmlOutput("kegg_effect")
                             ),
                             tabPanel("KEGG test results",
                                      htmlOutput("kegg_test")
                             ),
                             tabPanel("COG barplots",
                                      htmlOutput("cog_barplot_agg"),
                                      htmlOutput("cog_barplot_ind")
                             ),
                             tabPanel("COG boxplots",
                                      htmlOutput("cog_boxplot_agg"),
                                      htmlOutput("cog_boxplot_ind")
                             ),
                             tabPanel("COG effect size",
                                      htmlOutput("cog_effect")
                             ),
                             tabPanel("COG test results",
                                      htmlOutput("cog_test")
                             )
                      )
                  )
                )
        ),
        tabItem(tabName="subtype_discovery",
                fluidRow(
                  box(width=3,
                      h2("Subtype discovery"),
                      selectInput("subtype_method", "Subtype discovery method", choices=c("PAM", "DMM"), selected="PAM"),
                      selectInput("subtype_distance", "Distance used:", choices=c("UniFrac", "WUniFrac", "GUniFrac", "BC", "JS"), selected="UniFrac"),
                      checkboxGroupInput("assessment", "Assessment statistics:", choices=c("Gap", "ASW"), selected="Gap"),
                      actionButton("run_subtype", label="Submit"),
                      verbatimTextOutput("subtype_text")
                  ),
                  box(width=9,
                      tabBox(width=9,
                             tabPanel("Gap statistic",
                                      htmlOutput("gap_statistic")
                             ),
                             tabPanel("Avg. silhouette width",
                                      htmlOutput("silhouette_width")
                             ),
                             tabPanel("PCoA on UniFrac",
                                      htmlOutput("pcoa_unifrac")
                             ),
                             tabPanel("Cluster-specific taxa barplot",
                                      htmlOutput("cluster_barplot")
                             ),
                             tabPanel("Cluster-specific taxa boxplot",
                                      htmlOutput("cluster_boxplot")
                             ),
                             tabPanel("Cluster-specific effect size",
                                      htmlOutput("cluster_effect")
                             ),
                             tabPanel("Association test",
                                      htmlOutput("cluster_association")
                             )
                      )
                  )
                )
        ),
        tabItem(tabName="network_analysis",
                fluidRow(
                  box(width=3,
                      h2("Parameters"),
                      selectInput("graph_layout", "Graph layout:", c("Automatic" = "layout_nicely", "Bipartite" = "layout.bipartite",
                                                                     "Fruchterman-Reingold" = "layout_with_fr","Kamada-Kawai" = "layout_with_kk"), 
                                  selected="layout_nicely")
                  ),
                  box(width=9,
                      tabBox(width=NULL,
                             tabPanel("SpiecEasi output",
                                      plotOutput("spiec_easi", width="100%")
                             )
                      )
                  )
                )
        ),
        tabItem(tabName="Report",
                fluidRow(
                  box(width=3,
                      h2("Parameters"),
                      actionButton("run_report", label="Generate report"),
                      verbatimTextOutput("report_text"),
                      downloadLink("downloadData", "Download")
                  )
                ))
      )
    )
    
  ),
  
  server = function(input, output, session){
    
    data <- reactiveValues(val=NULL)
    data.rff <- reactiveValues(val=NULL)
    dist <- reactiveValues(val=NULL)
    dist.rff <- reactiveValues(val=NULL)
    phylo <- reactiveValues(val=NULL)
    phylo.rff <- reactiveValues(val=NULL)
    
    tables <- reactiveValues(phy.prev=NULL, phy.abund = NULL, fam.prev = NULL, fam.abund = NULL, gen.prev = NULL, gen.abund = NULL)
    
    alpha_results <- reactiveValues(rarefy_curve=NULL,boxplot=NULL,stats=NULL)
    
    beta <- reactiveValues(ord=NULL,clust=NULL,barplot=NULL,boxplot=NULL,permanova=NULL,mirkat=NULL,disper=NULL)
    
    diff_vis <- reactiveValues(val=NULL)
    
    samples_removed_vector <- reactiveValues(val=NULL)
    samples_kept_vector <- reactiveValues(val=NULL)
    
    source("ShinyStats.R")
    
    output$columns = renderUI({
      mydata = get(input$dataset)
      selectInput('columns2', 'Columns', names(mydata))
    })
    
    mapping_file <- reactive({input$mapping_file})
    df <- reactive({read.csv(mapping_file()$datapath, sep="\t")})
    
    
    #Selection
    output$maptest = renderUI({
      if(is.null(mapping_file())){
        helpText("Upload mapping file to use this functionality")
      }else{
        selectInput("filter_var", 'Filter by variable', c(names(df()), "Pick one"), "Pick one")
      }
    })
    output$maptest2 = renderUI({
      if(is.null(input$filter_var) || input$filter_var == "Pick one"){
      }else{
        dat <- df()
        cate <- input$filter_var
        DF <- data.frame(v1 = c(1,2,3,2), v2 = c("a","a","b","b"))
        checkboxGroupInput("filter_var_categories", 'Select categories to include in analysis. Unselected data will be filtered out.', as.character(unique(dat[[cate]])))
      }
    })
    output$maptest3 = renderUI({
      if(is.null(mapping_file())){
        helpText("Upload mapping file to use this functionality")
      }else{
        selectInput("sample_name", 'Sample name category (optional, only if using textbox below)', c(names(df()), "Pick one"), "Pick one")
      }
    })
    

    output$maptest6 = renderUI({
      if(is.null(input$covariate) || input$covariate == "Pick one"){
      }else{
        dat <- df()
        cate <- input$covariate
        DF <- data.frame(v1 = c(1,2,3,2), v2 = c("a","a","b","b"))
        checkboxGroupInput("selected_covariate", 'Pick category corresponding to covariates', as.character(unique(dat[[cate]])))
      }
    })
    
    output$output_text <- reactive({
      if(is.null(input$sample_name) || input$sample_name == "Pick one"){
        "Upload mapping file and select sample name category"
      }else{
        remove <- unlist(strsplit(input$input_text, split="\n"))
        sam <- input$sample_name
        dat <- df()
        names <- dat[[sam]]
        
        setdiff(remove, dat[[sam]])
      }
    })
    
    output$select_var = renderUI({
      if(is.null(mapping_file())){
        helpText("Upload mapping file to use this functionality")
      }else{
        selectInput("category", 'Variable of interest:', c("Pick one" = "", names(df())), "Pick one", selectize = TRUE)
      }
    })
    
    output$select_covars = renderUI({
      if(is.null(mapping_file())){
        helpText("Upload mapping file to use this functionality")
      }else{
        selectizeInput("covars", "Select covariates:", choices = c(names(df())), multiple=TRUE)
      }
    })
    
    output$select_subject = renderUI({
      if(is.null(mapping_file())){
        helpText("Upload mapping file to use this functionality")
      }else{
        selectInput("subject_var", 'Subject variable:', c("Pick one" = "", names(df())), "Pick one", selectize = TRUE)
      }
    })
    
    output$select_cat_vars = renderUI({
      if(is.null(mapping_file())){
        helpText("Upload mapping file to use this functionality")
      }else{
        selectizeInput("selected_cat_vars", "Select categorical variables:", choices = c(names(df())), multiple=TRUE)
      }
    })
    
    output$select_ord_vars = renderUI({
      if(is.null(mapping_file())){
        #helpText("Upload mapping file to use this functionality")
      }else{
        selectizeInput("selected_ord_vars", "Select ordinal variables:", choices = c(names(df())), multiple=TRUE)
      }
    })
    
    output$select_con_vars = renderUI({
      if(is.null(mapping_file())){
        #helpText("Upload mapping file to use this functionality")
      }else{
        selectizeInput("selected_con_vars", "Select continuous variables:", choices = c(names(df())), multiple=TRUE)
      }
    })
    
    #Alpha diversity
    output$alpha_measures = renderUI({
      measures = c("Observed", "Chao1", "Shannon", "InvSimpson")
      checkboxGroupInput("a_measures", "Alpha diversity measures", choices = measures, selected = measures)
    })
    
    #Beta diversity
    output$beta_measures = renderUI({
      measures = c("BC", "GUniFrac", "UniFrac", "WUniFrac")
      checkboxGroupInput("b_measures", "Beta diversity measures", choices = measures, selected = measures)
    })
    
    #Taxa diversity
    output$vis_level = renderUI({
      cols = c("Phylum","Class", "Order", "Family", "Genus", "Species")
      checkboxGroupInput("vis_level", "Visualization levels:", choices = cols, selected = cols)
    })
    output$func_vis_level = renderUI({
      cols = c("Phylum","Class", "Order", "Family", "Genus", "Species")
      checkboxGroupInput("vis_level", "Visualization levels:", choices = cols, selected = cols)
    })

    #submit button
    submit <- eventReactive(input$create_dataset,{
      num.var = input$selected_con_vars
      selection = input$r_select
      if(selection == ""){
        selection = NULL
      }
      n <- 8
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Creating dataset...", value = 0)
      progress$inc(1/n, detail = paste("Loading packages..."))
      progress$inc(1/n, detail = paste("Loading data..."))
      print(num.var)
      
      mappingFile <- input$mapping_file
      biomFile <- input$biom_file
      treeFile <- input$tree_file
      koFile <- input$kegg_file
      cogFile <- input$cog_file
      
      koAnn <- file.path(getwd(), "kegg.map.RData")
      data.obj <- load_data(otu.file=biomFile$datapath, map.file=mappingFile$datapath, tree.file=treeFile$datapath, ko.file=koFile$datapath, cog.file=cogFile$datapath, ko.ann.file=koAnn, meta.sep='\t', filter.no = input$otu_filter_count)
      data.obj.raw <- data.obj
      print(sum(colSums(data.obj$otu.tab)))
      print(sum(colSums(data.obj$abund.list[[1]])))
      progress$inc(1/n, detail = paste("Removing bad samples and OTUs..."))
      high_depth_samples <- rownames(data.obj$meta.dat)[colSums(data.obj$otu.tab) >= input$filter_dep]
      low_depth_samples <- rownames(data.obj$meta.dat)[colSums(data.obj$otu.tab) < input$filter_dep]
      usr_excl_samples <- unlist(strsplit(input$input_text, split="\n")) #from user-specified sample names
      bad_samples <- c(low_depth_samples, usr_excl_samples)
      samples_removed_vector$val <- colSums(data.obj$otu.tab)[low_depth_samples]
      good_samples <- rownames(data.obj$meta.dat)[! rownames(data.obj$meta.dat) %in% bad_samples ]
      samples_kept_vector$val <- colSums(data.obj$otu.tab)[good_samples]
      data.obj <- subset_data(data.obj, good_samples)
      print(sum(colSums(data.obj$otu.tab)))
      print(sum(colSums(data.obj$abund.list[[1]])))
      #low_freq_otus <- rownames(data.obj$otu.tab)[rowSums(data.obj$otu.tab != 0)  <= input$otu_filter_freq]
      #print(low_freq_otus)
      low_count_otus <- rownames(data.obj$otu.tab)[rowSums(data.obj$otu.tab) <= input$otu_filter_count]
      low_count_raw <- rownames(data.obj.raw$otu.tab)[rowSums(data.obj.raw$otu.tab) <= input$otu_filter_count]
      #print(low_count_otus)
      #bad_otus <- c(low_freq_otus, low_count_otus)
      #good_otus <- rownames(data.obj$otu.tab)[!rownames(data.obj$otu.tab) %in% bad_otus]
      good_otus <- rownames(data.obj$otu.tab)[!rownames(data.obj$otu.tab) %in% low_count_otus]
      good_otus_raw <- rownames(data.obj.raw$otu.tab)[!rownames(data.obj.raw$otu.tab) %in% low_count_raw]
      data.obj$otu.tab <- data.obj$otu.tab[good_otus,]
      data.obj.raw$otu.tab <- data.obj.raw$otu.tab[good_otus_raw,]
      print(sum(colSums(data.obj$otu.tab)))
      print(sum(colSums(data.obj$abund.list[[1]])))
      print(paste0(length(low_count_otus), " OTUs and ", length(bad_samples), " samples removed"))
      
      dist.obj <- construct_distance(data.obj)
      dist.obj.raw <- construct_distance(data.obj.raw)
      data.obj.rff <- load_data(otu.file=biomFile$datapath, map.file=mappingFile$datapath, tree.file=treeFile$datapath, rff=TRUE, dep=input$rare_dep, ko.file=koFile$datapath, cog.file=cogFile$datapath, ko.ann.file=koAnn,meta.sep='\t')
      dist.obj.rff <- construct_distance(data.obj.rff)
      
      save(data.obj, data.obj.rff, dist.obj, dist.obj.rff, file='Data.RData')
      
      progress$inc(1/n, detail = paste("Subsetting data..."))
      
      filter_var = input$filter_var
      print(filter_var)
      filter_var_cat = input$filter_var_categories
      print(filter_var_cat)
      print(selection)
      if (is.null(filter_var_cat)) {
        if (is.null(selection)) {
          filIDs <- rownames(data.obj$meta.dat)
          filIDs.rff <- rownames(data.obj.rff$meta.dat)
        }else{
          filIDs <- rownames(data.obj$meta.dat)[eval(parse(text=selection), envir=data.obj$meta.dat)]
          filIDs.rff <- rownames(data.obj.rff$meta.dat)[eval(parse(text=selection), envir=data.obj.rff$meta.dat)]
        }
        # remove NA	
        filIDs <- intersect(filIDs,  rownames(data.obj$meta.dat)[!is.na(data.obj$meta.dat)])
        filIDs.rff <- intersect(filIDs.rff,  rownames(data.obj.rff$meta.dat)[!is.na(data.obj.rff$meta.dat)])
      }else{
        if(is.null(selection)){
          filIDs <- rownames(data.obj$meta.dat)[data.obj$meta.dat[, filter_var] %in% filter_var_cat]
          filIDs.rff <- rownames(data.obj.rff$meta.dat)[data.obj.rff$meta.dat[, filter_var] %in% filter_var_cat]
        }else{
          print(paste0(is.null(selection)))
          print(paste0(selection))
          filIDs <- rownames(data.obj$meta.dat)[data.obj$meta.dat[, filter_var] %in% filter_var_cat & eval(parse(text=selection), envir=data.obj$meta.dat)] 
          filIDs.rff <- rownames(data.obj.rff$meta.dat)[data.obj.rff$meta.dat[, filter_var] %in% filter_var_cat & eval(parse(text=selection), envir=data.obj.rff$meta.dat)]
        }
        # remove NA	
        filIDs <- intersect(filIDs,  rownames(data.obj$meta.dat)[!is.na(data.obj$meta.dat[, filter_var])])
        filIDs.rff <- intersect(filIDs.rff,  rownames(data.obj.rff$meta.dat)[!is.na(data.obj.rff$meta.dat[, filter_var])])
      }
      print(sum(colSums(data.obj$otu.tab)))
      print(sum(colSums(data.obj$abund.list[[1]])))
      data.obj <- subset_data(data.obj, filIDs)
      data.obj.rff <- subset_data(data.obj.rff, filIDs.rff)
      dist.obj <- subset_dist(dist.obj, filIDs)
      dist.obj.rff <- subset_dist(dist.obj.rff, filIDs.rff)
      print(sum(colSums(data.obj$otu.tab)))
      print(sum(colSums(data.obj$abund.list[[1]])))
      #NOTE: Find somewhere to put this code later -- this will be used for covariates / variable types
      #if (grp.name %in% num.var) {
      #  data.obj$meta.dat[, grp.name] <- as.numeric(data.obj$meta.dat[, grp.name])
      #  data.obj.rff$meta.dat[, grp.name] <- as.numeric(data.obj.rff$meta.dat[, grp.name])
      #} else {
      #  data.obj$meta.dat[, grp.name] <- factor(data.obj$meta.dat[, grp.name], levels=grp.level.use)
      #  data.obj.rff$meta.dat[, grp.name] <- factor(data.obj.rff$meta.dat[, grp.name], levels=grp.level.use)
      #}
      
      colnames(data.obj$meta.dat) <- gsub("^\\s+|\\s+$", "", colnames(data.obj$meta.dat))
      colnames(data.obj.rff$meta.dat) <- gsub("^\\s+|\\s+$", "", colnames(data.obj.rff$meta.dat))
      print(colSums(data.obj$otu.tab))
      print(colSums(data.obj$abund.list[[1]]))
      progress$inc(1/n, detail = paste("Creating phyloseq object..."))
      
      phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
                            tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
      phylo.obj.rff <- phyloseq(otu_table(data.obj.rff$otu.tab, taxa_are_rows=T), phy_tree(data.obj.rff$tree), 
                                tax_table(data.obj.rff$otu.name), sample_data(data.obj.rff$meta.dat))
      data$val = data.obj
      data.rff$val = data.obj.rff
      dist$val = dist.obj
      dist.rff$val = dist.obj.rff
      phylo$val = phylo.obj
      phylo.rff$val = phylo.obj
      
      progress$inc(1/n, detail = paste("Dataset created."))
      save(data.obj.raw, data.obj, data.obj.rff, dist.obj.raw, dist.obj, dist.obj.rff, file='Data.RData')
      save(phylo.obj, file='Phylo.RData')
      save(phylo.obj.rff, file='PhyloRar.RData')
      print("Dataset created!")
    })
    
    output$status <- renderText({
      submit()
    })
    
    observeEvent(input$category,{
      dat <- df()
      grp.level.use <- unique(dat[[input$category]])
      if (input$category %in% input$selected_con_vars) {
        data$val$meta.dat[, input$category] <- as.numeric(data$val$meta.dat[, input$category])
        data.rff$val$meta.dat[, input$category] <- as.numeric(data.rff$val$meta.dat[, input$category])
      } else {
        data$val$meta.dat[, input$category] <- factor(data$val$meta.dat[, input$category], levels=grp.level.use)
        data.rff$val$meta.dat[, input$category] <- factor(data.rff$val$meta.dat[, input$category], levels=grp.level.use)
      }
    })
    
    
    submit_summary <- eventReactive(input$summary_stats,{
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      n <- 3
      progress$set(message = "Run summary statistics...", value = 0)
      progress$inc(1/n, detail = paste("Performing summary analysis..."))
      otu.tab <- data$val$otu.tab
      prev = input$prev / 100
      abund = input$abund / 100

      otu.abund <- rowSums(otu.tab)
      sam.abund <- colSums(otu.tab)
      otu.prev <- rowSums(otu.tab!=0)/ncol(otu.tab)
      
      otu.abund <- otu.abund[otu.abund >= 1]
      sam.abund <- sam.abund[sam.abund >= 1]
      
      phy.abund <- data$val$abund.list[['Phylum']]
      fam.abund <- data$val$abund.list[['Family']]
      gen.abund <- data$val$abund.list[['Genus']]
      
      phy.prev <- rowSums(phy.abund != 0) / ncol(phy.abund)
      fam.prev <- rowSums(fam.abund != 0) / ncol(phy.abund)
      gen.prev <- rowSums(gen.abund != 0) / ncol(phy.abund)
      
      phy.abund <- rowMeans(t(t(phy.abund) / sam.abund))
      fam.abund <- rowMeans(t(t(fam.abund) / sam.abund))
      gen.abund <- rowMeans(t(t(gen.abund) / sam.abund))
      
      phy.prev <- sort(phy.prev, decr=T)
      fam.prev <- sort(fam.prev, decr=T)
      gen.prev <- sort(gen.prev, decr=T)
      
      phy.abund <- sort(phy.abund, decr=T)
      fam.abund <- sort(fam.abund, decr=T)
      gen.abund <- sort(gen.abund, decr=T)
      
      tables$phy.prev <- round(phy.prev[phy.prev >= 0.05] * 100, 2)
      tables$fam.prev <- round(fam.prev[fam.prev >= 0.05] * 100, 2)
      tables$gen.prev <- round(gen.prev[gen.prev >= 0.05] * 100, 2)
      tables$phy.abund <- round(phy.abund[phy.abund >= 0.05] * 100, 2)
      tables$fam.abund <- round(fam.abund[fam.abund >= 0.05] * 100, 2)
      tables$gen.abund <- round(gen.abund[gen.abund >= 0.05] * 100, 2)
      
      progress$inc(1/n, detail = paste("Generating heatmaps..."))
      progress$inc(1/n, detail = paste("Generating barplots..."))
      print("Done!")
    })
    
    observeEvent(input$summary_stats,{
      output$cov_dist <- renderPlot({
        sam.abund <- colSums(data$val$otu.tab)
        sam.abund <- sam.abund[sam.abund >= 1]
        ggplot2::ggplot(data=data.frame(x=sam.abund), aes(x=x)) + geom_histogram(col='black', fill='gray')  + ylab('Frequency') + xlab('Sequencing depth') + theme_bw()
      })
      output$cov_boxplot <- renderPlot({
        otu.tab <- data$val$otu.tab
        map <- data$val$meta.dat
        colnames(otu.tab) <-  map[[input$category]]
        df <- data.frame(Group=names(colSums(otu.tab)), coverage=colSums(otu.tab))
        ggplot2::ggplot(df, aes(x=Group, y=log10(coverage), col=Group)) + geom_boxplot(position=position_dodge(width=0.75), outlier.colour = NA) +
          geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) + theme_bw()
      })
      
      output$phy.prev <- DT::renderDataTable({
        reshape2::melt(tables$phy.prev, value.name="Prevalence (%)")
      }, caption="Prevalent Phyla")
      output$fam.prev <- DT::renderDataTable({
        reshape2::melt(tables$fam.prev, value.name="Prevalence (%)")
      }, caption="Prevalent Families")
      output$gen.prev <- DT::renderDataTable({
        reshape2::melt(tables$gen.prev, value.name="Prevalence (%)")
      }, caption="Prevalent Genus")
    
      output$phy.abund <- DT::renderDataTable({
        reshape2::melt(tables$phy.abund, value.name="Abundance (%)")
      }, caption="Abundant Phyla")
      output$fam.abund <- DT::renderDataTable({
        reshape2::melt(tables$fam.abund, value.name="Abundance (%)")
      }, caption="Abundant Families")
      output$gen.abund <- DT::renderDataTable({
        reshape2::melt(tables$gen.abund, value.name="Abundance (%)")
      }, caption="Abundant Genus")
    })
    
    observeEvent({
      input$summary_stats
      input$heatmap_type
    },{
      output$summary_heatmap <- renderIheatmap({
        prop <- prop.table(data$val$abund.list$Genus,2)
        if(input$heatmap_type == "Proportional"){
          col.scheme = c("white", brewer.pal(11, "Spectral"))
          minp <- min(prop[prop!=0]/1.1)
          prop[prop==0] <- minp
          prop <- log10(prop)
        }else if(input$heatmap_type == "Binary"){
          col.scheme <- c("lightyellow", "red")
          prop[, ] <- as.numeric(prop != 0)
        }else if(input$heatmap_type == "Ranked"){
          col.scheme <- c('white', colorRampPalette(c("green", "black", "red"))(ncol(prop)-1))
          prop <- t(apply(prop, 1, function(x) {
            temp <- rank(x[x!=0])
            s <- (ncol(prop) - 1) / (max(temp) - min(temp))
            temp <- 1 + (temp - min(temp)) * s
            x[x!=0] <- temp
            x
          }))
        }
        phy <- sapply(strsplit(rownames(prop), ";"), function(x) x[1])
        main_heatmap(prop, colors=col.scheme) %>% 
          add_row_annotation(phy) %>%
          add_col_annotation(as.data.frame(data$val$meta.dat[,input$category, drop=FALSE])) %>%
          add_row_clustering() %>%
          add_col_clustering()
      })
    })
    
    
    
    observeEvent({
      input$summary_stats
      input$barplot_level
    },{
      prop <- prop.table(data$val$abund.list[[input$barplot_level]],2)
      prop.m <- melt(prop[rev(order(rowMeans(prop))),])
      prop.m$factor1 <- data$val$meta.dat[match(prop.m$Var2, rownames(data$val$meta.dat)), input$category]
      
      output$summary_barplot1 <- renderPlot({
        ggplot(prop.m, aes(factor1, value, fill = Var1, key=Var1) ) +
        geom_bar(stat="identity") +
        guides(fill=FALSE) + scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(prop.m$Var1))))
      })
      
      output$summary_barplot2 <- renderPlot({
        ggplot(prop.m, aes(Var2, value, fill = Var1) ) +
          geom_bar(stat="identity") +
          guides(fill=FALSE) + facet_grid(~factor1, scales="free", space="free_x") + scale_fill_manual(values = colorRampPalette(brewer.pal(12, "Set3"))(length(unique(prop.m$Var1))))
      })
    }, ignoreInit = TRUE)
    
    output$summary_status <- renderText({
      submit_summary()
    })
   
    submit_alpha <- eventReactive(input$run_alpha,{
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Making plot", value = 0)
      n <- 3
      
      progress$inc(1/n, detail = paste("Generating rarefication curve...", 1))
      print(input$category)
      alpha_results$rarefy_curve <- generate_rarefy_curve(data$val, phylo$val, grp.name=input$category, depth=input$alpha_rare_dep, iter.no=input$rare_iter)
      progress$inc(1/n, detail = paste("Generating alpha diversity boxplots...", 2))
      alpha_results$boxplot <- generate_alpha_boxplot(data$val, phylo$val, depth=input$alpha_rare_dep, grp.name=input$category, strata=NULL)
      progress$inc(1/n, detail = paste("Perfoming association tests...", 3))
      alpha_results$stats <- perform_alpha_test2(data$val, depth=input$alpha_rare_dep, iter.no=input$rare_iter, grp.name=input$category, adj.name=NULL)
      print("Done!")
    })
    
    output$alpha_text <- renderText({
      submit_alpha()
    })
    
    observeEvent(input$run_alpha,{
      output$rarefy_curve <- renderPlot({
        alpha_results$rarefy_curve
      })
      output$rarefy_boxplot <- renderPlot({
        alpha_results$boxplot
      })
      
      output$alpha_association_tab1 <- renderUI({
        HTML(print(xtable(summary(alpha_results$stats$fitted.obj$Observed)$coefficients, caption="Observed test results:"), 
                   type="html",
                   caption.placement="top",
                   html.table.attributes='class="data table table-bordered table-condensed"'))
      })
      output$alpha_association_tab2 <- renderUI({
        HTML(print(xtable(anova(alpha_results$stats$fitted.obj$Observed), caption="Observed ANOVA results:"), 
                   type="html",
                   caption.placement="top",
                   html.table.attributes='class="data table table-bordered table-condensed"'))
      })
      output$alpha_association_tab3 <- renderUI({
        HTML(print(xtable(summary(alpha_results$stats$fitted.obj$Chao1)$coefficients, caption="Chao1 test results:"), 
                   type="html",
                   caption.placement="top",
                   html.table.attributes='class="data table table-bordered table-condensed"'))
      })
      output$alpha_association_tab4 <- renderUI({
        HTML(print(xtable(anova(alpha_results$stats$fitted.obj$Chao1),caption="Chao1 ANOVA results:"), 
                   type="html",
                   caption.placement="top",
                   html.table.attributes='class="data table table-bordered table-condensed"'))
      })
      output$alpha_association_tab5 <- renderUI({
        HTML(print(xtable(summary(alpha_results$stats$fitted.obj$Shannon)$coefficients,caption="Shannon test results:"), 
                   type="html",
                   caption.placement="top",
                   html.table.attributes='class="data table table-bordered table-condensed"'))
      })
      output$alpha_association_tab6 <- renderUI({
        HTML(print(xtable(anova(alpha_results$stats$fitted.obj$Shannon),caption="Shannon ANOVA results:"), 
                   type="html",
                   caption.placement="top",
                   html.table.attributes='class="data table table-bordered table-condensed"'))
      })
      output$alpha_association_tab7 <- renderUI({
        HTML(print(xtable(summary(alpha_results$stats$fitted.obj$InvSimpson)$coefficients,caption="InvSimpson test results:"), 
                   type="html",
                   caption.placement="top",
                   html.table.attributes='class="data table table-bordered table-condensed"'))
      })
      output$alpha_association_tab8 <- renderUI({
        HTML(print(xtable(anova(alpha_results$stats$fitted.obj$InvSimpson), caption="InvSimpson ANOVA results:"), 
                   type="html",
                   caption.placement="top",
                   html.table.attributes='class="data table table-bordered table-condensed"'))
      })
    })
    
    submit_beta <- eventReactive(input$run_beta,{
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin beta diversity analysis", value = 0)
      n <- 7
      progress$inc(1/n, detail = paste("Generating ordination plots..."))

      beta$ord <- generate_ordination2(data.rff$val, dist.rff$val, grp.name=input$category, strata=NULL, dist.names = input$b_measures)
      
      progress$inc(1/n, detail = paste("Generating boxplots..."))

      beta$boxplot <- generate_distance_boxplot(data.rff$val, dist.rff$val, grp.name=input$category, within=T, strata=NULL)
      
      progress$inc(1/n, detail = paste("Performing permanova test..."))
      beta$permanova <- perform_permanova_test(data.rff$val, dist.rff$val, grp.name=input$category, adj.name=NULL, strata=NULL)
      
      progress$inc(1/n, detail = paste("Performing mirkat test..."))

      beta$mirkat <- perform_mirkat_test(data.rff$val, dist.rff$val, grp.name=input$category, adj.name=NULL)    # Could not handle correlation
      
      progress$inc(1/n, detail = paste("Performing beta dispersion test..."))
      beta$disper <- perform_betadisper_test(data.rff$val, dist.rff$val, grp.name=input$category)
      
      print("Done!")
    })
    
    output$beta_text <- renderText({
      submit_beta()
    })
    
    observeEvent(input$run_beta,{
      output$ordination <- renderPlotly({
        ggplotly(beta$ord)
      })
      output$distance_comparison_boxplot <- renderPlot({
        beta$boxplot
      })
      output$permanova <- renderUI({
        out <- permanova_tab()
        div(HTML(as.character(out)),class="shiny-html-output")
      })
      output$mirkat <- renderUI({
        out <- mirkat_tab()
        div(HTML(as.character(out)),class="shiny-html-output")
      })
      output$betadisper <- renderUI({
        out <- betadisper_tab()
        div(HTML(as.character(out)),class="shiny-html-output")
      })
    })
    
    permanova_tab <- function(){
      measures <- input$b_measures
      tables <- list()
      tables <- lapply(measures, function(x){
        print(xtable(beta$permanova$permanova.obj[[x]], caption=paste(x, "Permanova")), 
              type="html", 
              html.table.attributes='class="data table table-bordered table-condensed"', 
              caption.placement="top")
      })
      G_caption <- paste('PERMANOVA G test combining ', paste(measures, collapse=','))
      tables[[as.character(length(measures)+1)]] <- print(xtable(beta$permanova$permanovaG.obj, caption=G_caption),
                                                          type="html", 
                                                          html.table.attributes='class="data table table-bordered table-condensed"', 
                                                          caption.placement="top")
      all <- lapply(tables, paste0)
      return(all)
    }
    
    mirkat_tab <- function(){
      mirkat <- cbind(beta$mirkat$indiv, beta$mirkat$omni)
      colnames(mirkat) =  c("UniFrac", "GUniFrac", "WUniFrac", "BC", "Omnibus")
      all <- print(xtable(mirkat, caption=paste("P-values for MiRKAT test combining UniFrac,GUniFrac,WUniFrac,BC"), digits=4),
                   type="html", 
                   html.table.attributes='class="data table table-bordered table-condensed"', 
                   caption.placement="top")
      return(all)
    }
    
    betadisper_tab <- function(){
      measures <- input$b_measures
      tables <- list()
      tables <- lapply(measures, function(x){
        print(xtable(beta$disper[[x]], caption=paste(x, "BETADISPER")), 
              type="html", 
              html.table.attributes='class="data table table-bordered table-condensed"', 
              caption.placement="top")
      })
      all <- lapply(tables, paste0)
      return(all)
    }
    
    output$taxa_text <- renderText({
      submit_taxa()
      print("Done!")
    })
    
    submit_taxa <- eventReactive(input$run_taxa,{
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin taxa diversity analysis", value = 0)
      n <- 4
      progress$inc(1/n, detail = paste("Performing differential taxa analysis..."))
      set.seed(123)
      diff.obj.rff <- perform_differential_analysis(data.rff$val, 
                                                    grp.name=input$category, 
                                                    adj.name=NULL, 
                                                    taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species'),
                                                    method=input$taxa_method, 
                                                    mt.method=input$mult_test, 
                                                    subject=NULL,
                                                    cutoff=input$sig_level / 100, 
                                                    prev=input$taxa_prev / 100, 
                                                    minp=input$taxa_abund / 100, 
                                                    ann=input$taxa_method)
      progress$inc(1/n, detail = paste("Creating visualizations for differential analysis..."))
      diff_vis$val <- visualize_differential_analysis(data.rff$val, 
                                      diff.obj.rff, 
                                      grp.name=input$category, 
                                      taxa.levels=input$vis_level, 
                                      mt.method=input$mult_test, 
                                      cutoff=input$sig_level / 100, 
                                      ann=input$taxa_method)
      progress$inc(1/n, detail = paste("Beginning LEfSe analysis. Creating format..."))
      #create_lefse_format(data.rff$val, 
      #                    diff.obj.rff, 
      #                    grp.name=input$category, 
      #                    cutoff=input$sig_level / 100, 
      #                    prev=input$taxa_prev / 100, 
      #                    minp=input$taxa_abund / 100, 
      #                    mt.method=input$mult_test)
      progress$inc(1/n, detail = paste("Performing LEfSe analysis..."))
      #system("source runENV.sh")
      save(diff.obj.rff, file="DiffData.RData")
      #perform_lefse_analysis(data.rff$val, grp.name=input$category, ann=input$category)
    })
    
    observeEvent(input$run_taxa,{
      
      output$taxa_boxplots <- renderPlot({
        diff_vis$val$boxplot_aggregate
      })
      output$taxa_barplots <- renderPlot({
        diff_vis$val$barplot_aggregate
      })
      output$effect_size <- renderPlot({
        diff_vis$val$effect_size
      })
      output$taxa_biplot <- renderPlot({
        diff_vis$val$biplot   
      })
      output$taxa_prop_heatmap <- renderIheatmap({
        diff_vis$val$prop_heatmap
      })
      output$taxa_rank_heatmap <- renderIheatmap({
        diff_vis$val$rank_heatmap
      })
      #output$cladogram <- renderText({
      #  name <- paste0('<iframe style="height:600px; width:900px" src="plots/LefSe_', input$category, '/cladogram.pdf"></iframe>')
      #  return(name)
      #})
      output$taxa_test_results <- renderUI({
        out <- taxa_test_tab()
        div(HTML(as.character(out)),class="shiny-html-output")
      })
    })
    
    taxa_test_tab <- function(){
      mytab <- read.csv(paste0("Taxa_DifferentialAbundanceAnalysis_AllLevels_", input$mult_test, "_", input$sig_level/100, "_", input$taxa_method,".csv"))
      names(mytab)[1] <- "Taxa"
      table <- print(xtable(mytab[,c(1,2,3,4,5,8)], caption="Differential abundance analysis for all levels"),
                     type="html", 
                     html.table.attributes='class="data table table-bordered table-condensed"', 
                     caption.placement="top")
      return(table)
    }
    
    output$pred_text <- renderText({
      submit_pred()
      print("Done!")
    })
    
    submit_pred <- eventReactive(input$run_pred,{
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Running prediction modeling", value = 1/2)
      predictionRF(data$val,  resp.name=input$category, taxa.level='Genus', boruta.leve=input$boruta_level, B=input$bootstrap_num, prev=input$pred_prev / 100, minp=input$pred_abund / 100)
    })
    
    observeEvent(input$run_pred,{
      output$classification_error <- renderImage({
        input$run_pred
        list(src = "Taxa_Random_forest_misclassification_barplot_Genus_.png",
             contentType = 'image/png',
             alt = "Misclassification barplot")
      }, deleteFile = FALSE)
      
      output$bootstrap_roc_genus <- renderImage({
        input$run_pred
        list(src = "Taxa_Random_forest_ROC_Genus_.png",
             contentType = 'image/png',
             alt = "Bootstrap validated ROC curve at Genus level")
      }, deleteFile = FALSE)
      
      output$boostrap_roc_species <- renderImage({
        input$run_pred
        list(src = "Taxa_Random_forest_ROC_Species_.png",
             contentType = 'image/png',
             alt = "Bootstrap validated ROC curve at Species level")
      }, deleteFile = FALSE)
      
      output$boruta_features <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_Random_forest_Boruta_Feature_Selection_Genus_.pdf"></iframe>')
        return(name)
      })
      
      output$boruta_barplots_agg <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_Barplot_Aggregate_Genus_sqrt_BorutaFeatures_Tentative__.pdf"></iframe>')
        return(name)
      })
      
      output$boruta_barplots_ind <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_Barplot_Genus_P_BorutaFeatures_Tentative__.pdf"></iframe>')
        return(name)
      })
      
      output$boruta_boxplots <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_Boxplot_Genus_P_BorutaFeatures_Tentative__.pdf"></iframe>')
        return(name)
      })
      
      output$boruta_roc <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/BorutaFeatures_Tentative_ROC_Genus_0.632+.pdf"></iframe>')
        return(name)
      })
    })
    
    output$func_text <- renderText({
      submit_func()
      print("Done!")
    })
    
    submit_func <- eventReactive(input$run_func,{
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin functional diversity analysis", value = 0)
      n <- 8
      progress$inc(1/n, detail = paste("Generating KEGG boxplots..."))
      generate_taxa_boxplot(data.rff$val, grp.name=input$category, strata=NULL, rm.outlier=F, prev=input$func_prev / 100, minp=input$func_abund / 100, 
                            taxa.levels=c('KEGG_Metabolism'), ann='KEGG')
      progress$inc(1/n, detail = paste("Generating KEGG barplots..."))
      generate_taxa_barplot(data.rff$val, grp.name=input$category, strata=NULL, prev=input$func_prev / 100, minp=input$func_abund / 100,
                            taxa.levels=c('KEGG_Metabolism'), ann='KEGG')
      progress$inc(1/n, detail = paste("Performing KEGG differential analysis..."))
      diff.obj.rff <- perform_differential_analysis(data.rff$val, 
                                                    grp.name=input$category, 
                                                    adj.name=NULL,
                                                    method=input$func_method, 
                                                    mt.method=input$func_mult_test, 
                                                    prev=input$func_prev / 100, 
                                                    minp=input$func_abund / 100, 
                                                    cutoff=input$func_sig_level / 100, 
                                                    taxa.levels=c('KEGG_Metabolism'), 
                                                    ann=paste0('KEGG_', input$func_method))
      progress$inc(1/n, detail = paste("Generating KEGG visualizations..."))
      visualize_differential_analysis(data.rff$val, diff.obj.rff, grp.name=input$category, cutoff=input$func_sig_level / 100, taxa.levels=c('KEGG_Metabolism'), 
                                      ann='KEGG', scale='none', mt.method=input$func_mult_test)
      
      progress$inc(1/n, detail = paste("Generating COG boxplots..."))
      generate_taxa_boxplot(data.rff$val, grp.name=input$category, strata=NULL, rm.outlier=F, prev=input$func_prev / 100, minp=input$func_abund / 100,
                            taxa.levels=c('COG_Category2'), ann='COG')
      progress$inc(1/n, detail = paste("Generating COG barplots..."))
      generate_taxa_barplot(data.rff$val, grp.name=input$category, strata=NULL, prev=input$func_prev / 100, minp=input$func_abund / 100,
                            taxa.levels=c('COG_Category2'), ann='COG')
      progress$inc(1/n, detail = paste("Performing COG differential analysis..."))
      diff.obj.rff <- perform_differential_analysis(data.rff$val, 
                                                    grp.name=input$category, 
                                                    adj.name=NULL,
                                                    method=input$func_method, 
                                                    mt.method=input$func_mult_test, 
                                                    prev=input$func_prev / 100, 
                                                    minp=input$func_abund / 100,
                                                    cutoff=input$func_sig_level / 100,
                                                    taxa.levels=c('COG_Category2'), 
                                                    ann=paste0('COG_', input$func_method))
      progress$inc(1/n, detail = paste("Generating COG visualizations..."))
      visualize_differential_analysis(data.rff$val, diff.obj.rff, grp.name=input$category, cutoff=input$func_sig_level / 100, taxa.levels=c('COG_Category2'), 
                                      ann='COG', scale='none', mt.method=input$func_mult_test)
      ##save
    })
    observeEvent(input$run_func,{
      output$kegg_barplot_agg <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_AbundanceBarplot_none_fdr_0.1_KEGG_.pdf"></iframe>')
        return(name)
      })
      output$kegg_barplot_ind <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_Barplot_All_P_fdr_0.1_KEGG_.pdf"></iframe>')
        return(name)
      })
      output$kegg_boxplot_agg <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_AbundanceBoxplot_none_fdr_0.1_KEGG_.pdf"></iframe>')
        return(name)
      })
      output$kegg_boxplot_ind <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_Boxplot_All_P_fdr_0.1_KEGG_.pdf"></iframe>')
        return(name)
      })
      output$kegg_effect <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_logPBarplot_fdr_0.1_KEGG_.pdf"></iframe>')
        return(name)
      })
      output$kegg_test <- renderUI({
        out <- kegg_test_tab()
        div(HTML(as.character(out)),class="shiny-html-output")
      })
      output$cog_barplot_agg <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_AbundanceBarplot_none_fdr_0.1_COG_.pdf"></iframe>')
        return(name)
      })
      output$cog_barplot_ind <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_Barplot_All_P_fdr_0.1_COG_.pdf"></iframe>')
        return(name)
      })
      output$cog_boxplot_agg <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_AbundanceBoxplot_none_fdr_0.1_COG_.pdf"></iframe>')
        return(name)
      })
      output$cog_boxplot_ind <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_Boxplot_All_P_fdr_0.1_COG_.pdf"></iframe>')
        return(name)
      })
      output$cog_effect <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_logPBarplot_fdr_0.1_COG_.pdf"></iframe>')
        return(name)
      })
      output$cog_test <- renderUI({
        out <- cog_test_tab()
        div(HTML(as.character(out)),class="shiny-html-output")
      })
    })
    
    kegg_test_tab <- function(){
      mytab <- read.csv("Taxa_DifferentialAbundanceAnalysis_AllLevels_fdr_0.1_KEGG_perm.csv", head=TRUE, sep=",") 
      names(mytab)[1] <- "Taxa"
      table <- print(xtable(mytab[,c(1,2,3,4,5,8)], caption="Differential abundance analysis for all levels"),
                     type="html", 
                     html.table.attributes='class="data table table-bordered table-condensed"', 
                     caption.placement="top")
      return(table)
    }
    
    cog_test_tab <- function(){
      mytab <- read.csv("Taxa_DifferentialAbundanceAnalysis_AllLevels_fdr_0.1_COG_perm.csv", head=TRUE, sep=",") 
      names(mytab)[1] <- "Taxa"
      table <- print(xtable(mytab[,c(1,2,3,4,5,8)], caption="Differential abundance analysis for all levels"),
                     type="html", 
                     html.table.attributes='class="data table table-bordered table-condensed"', 
                     caption.placement="top")
      return(table)
    }
    
    output$subtype_text <- renderText({
      submit_subtype()
      print("Done!")
    })
    
    submit_subtype <- eventReactive(input$run_subtype,{
      progress <- shiny::Progress$new()
      on.exit(progress$close())
      progress$set(message = "Begin subtype analyis analysis. Beginning random forest predictions...", value = 0)
      n <- 2
      predictionRF(data$val,  resp.name=input$category, taxa.level='Genus', boruta.leve=input$boruta_level, B=input$bootstrap_num, prev=input$prev / 100, minp=input$abund / 100)
      progress$inc(1/n, detail = paste("Performing cluster analysis..."))
      perform_cluster_analysis(data$val, dist$val, dist.name='UniFrac', method='pam', stat='gap', 
                               grp.name=input$category, adj.name=NULL) 
    })
    
    observeEvent(input$run_subtype,{
      output$gap_statistic <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/cluster_assess_gap_statistic_UniFrac.pdf"></iframe>')
        return(name)
      })
      output$silhouette_width <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/cluster_assess_asw_statistic_UniFrac.pdf"></iframe>')
        return(name)
      })
      output$pcoa_unifrac <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Beta_diversity_ordination_cmd_ClusterUniFrac.pdf"></iframe>')
        return(name)
      })
      output$cluster_boxplot <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_AbundanceBoxplot_sqrt_fdr_0.01_Cluster_.pdf"></iframe>')
        return(name)
      })
      output$cluster_barplot <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_AbundanceBarplot_sqrt_fdr_0.01_Cluster_.pdf"></iframe>')
        return(name)
      })
      output$cluster_effect <- renderText({
        name <- paste0('<iframe style="height:600px; width:900px" src="plots/Taxa_DifferentialAbundance_logPBarplot_fdr_0.01_Cluster_.pdf"></iframe>')
        return(name)
      })
      output$cluster_association <- renderUI({
        out <- cluster_tab()
        div(HTML(as.character(out)),class="shiny-html-output")
      })
    })
    
    cluster_tab <- function(){
      lines <- readLines("Cluster_association_testUniFrac.txt")
      clusters <- grep("Test for enrichment", lines)
      fixed_colns <- character(5)
      start <- clusters + 1
      diff <- clusters[2] - start[1] - 1
      end <- start + diff
      
      tables <- list()
      
      for (i in 1:length(clusters)){
        tab <- read.table(text = lines[start[i]:end[i]], fill=TRUE, header=TRUE)
        colns <- strsplit (lines[start[i]], "\\s+")[[1]]
        fixed_colns <- c(" ", "Estimate", "Std. Error", "z value", "Pr(>|z|)")
        colnames(tab) <- fixed_colns
        cluster <- lines[clusters[i]]
        tables[[as.character(i)]] <- print(xtable(tab[,c(1,2,3,4,5)], caption=paste(cluster, "test results")),
                                           type="html", 
                                           html.table.attributes='class="data table table-bordered table-condensed"', 
                                           caption.placement="top")
      }
      all <- lapply(tables, paste0)
      return(all)
    }
    
    output$spiec_easi <- renderPlot({
      load("DiffData.RData")
      diff_genus <- gsub(".*;", "", rownames(diff.obj.rff$qv.list[['Genus']])[which(diff.obj.rff$qv.list[['Genus']] <= 0.1)])
      phylo_select <- subset_taxa(phylo.obj, Genus %in% diff_genus)
      
      categories <- unique(get_variable(phylo_select, VOI))
      
      networks <- lapply(categories, function(y){ 
        phylo_split <- prune_samples(get_variable(phylo_select, VOI)==y, phylo_select)
        se.mb <- spiec.easi(phylo_split, method='mb', lambda.min.ratio=1e-2, nlambda=20, icov.select.params=list(rep.num=50, verbose=TRUE, ncores=10))
        ig.mb <- adj2igraph(symBeta(getOptBeta(se.mb), mode='maxabs'), vertex.attr=list(name=tax_table(phylo_split)[,"Genus"]))
        set_vertex_attr(ig.mb, VOI, index=V(ig.mb1), y)
      })
      
      color_count <- length(diff_genus)
      getPalette = colorRampPalette(brewer.pal(length(diff_genus), "Set3"))
      
      disj <- disjoint_union(networks)
      V(disj2)$degree <- degree(disj)
      network_plot <- ggplot(ggnetwork(disj2, layout="fruchtermanreingold", by=VOI), aes(x = x, y = y, xend = xend, yend = yend)) + 
        geom_nodes(aes(color = vertex.names, size=degree)) + 
        geom_edges(aes(color = ifelse(weight > 0, 'green', 'red'))) + 
        facet_wrap(as.formula(paste("~", VOI))) + 
        theme_facet()
    })
    
    output$report_text <- renderText({
      submit_report()
      print("Done!")
    })
    
    submit_report <- eventReactive(input$run_report,{
      VOI <- input$category
      summary_table <- dplyr::select(data$val$meta.dat, VOI) %>% group_by_(VOI) %>% dplyr::summarize(n()) %>% knitr::kable()
      filter_dep <- input$filter_dep
      OTU_vector <- as.vector(rowSums(data$val$otu.tab))
      num_phyla <- unique(data$val$otu.name.full[,"Phylum"])
      num_family <- unique(data$val$otu.name.full[,"Family"])
      num_genus <- unique(data$val$otu.name.full[,"Genus"])
      samples_removed <- samples_removed_vector$val
      samples_kept <- samples_kept_vector$val
      perc_zero <- sum(colSums(data$val$otu.tab == 0))/(nrow(data$val$otu.tab)*ncol(data$val$otu.tab))*100
      num_rows <- length(unique(data$val$meta.dat[[VOI]]))
      rmarkdown::render("summary.Rmd", params = list(
        voi = VOI,
        table = summary_table,
        minreads = filter_dep,
        samples_removed = samples_removed,
        samples_kept = samples_kept,
        OTU_vector = OTU_vector,
        num_phyla = num_phyla,
        num_family = num_family,
        num_genus = num_genus,
        perc_zero = perc_zero
      ))
      rmarkdown::render("alpha_beta_diversity.Rmd", params = list(
        voi = VOI,
        num_rows = num_rows
      ))
      rmarkdown::render("taxa_diversity.Rmd", params = list(
        voi = VOI
      ))
      rmarkdown::render("predictive_modeling.Rmd", params = list(
        voi = VOI
      ))
      rmarkdown::render("functional_analysis.Rmd", params = list(
        voi = VOI
      ))
    })
    
    output$downloadData <- downloadHandler(
      filename = function() {
        paste("summary.html")
      },
      content <- function(file) {
        file.copy("summary.html", file)
      }
    )
  } 
)

