# MicrobiomeGPS

MicrobiomeGPS is an R/Shiny application which takes as input an OTU count table, phylogenetic tree and sample mapping file, with optional input for COG and KEGG tables. MicrobiomeGPS then performs numerous statistical tests and generates interactive visualizations and tables for different type of analysis such as exploratory data analysis, alpha- and beta- diversity analysis, taxonomic and functional differential abundance analysis, predictive modeling using random forests, community subtype analysis, and OTU network analysis. We include the most robust and powerful statistical methods developed recently for microbiome data. Specifically, we use linear model/mixed effects modeling for alpha-diversity, MiRKAT/PERMANOVA for beta-diversity, permutation-based FDR control for taxa/function data, tree-based random forests for predictive modeling, and SPIEC-EASI for OTU network analysis.  MicrobiomeGPS integrates covariate adjustment in each step and thus can address potential confounding effects due to technical, biological and clinical variables. It can also analyze repeated measurement data by taking into account the within-subject correlation structure. MicrobiomeGPS also allows users to import/export all parameters and download detailed HTML reports to facilitate collaboration and reproducibility. 

## Starting MicrobiomeGPS

### Via Github

Starting the MicrobiomeGPS Shiny app is quite simple. In an R session, simply enter:

```
shiny::runGitHub("SJohnsonMayo/MicrobiomeGPS")
```

###Via ShinyApps.io

MicrobiomeGPS is also currently deployed on a free shinyapps.io server: 

https://sjohnsonmayo.shinyapps.io/github_MicrobiomeGPS/


### Documentation

Watch this space for a manual and vignettes. 

### Issues and Updates

MicrobiomeGPS is still under active development. Upcoming features and any bugs can be found on the [issues page](https://github.com/SJohnsonMayo/MicrobiomeGPS/issues).
