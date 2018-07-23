load_data_new <- function (otu.file, map.file, tree.file=NULL, parseFunction=parse_taxonomy_greengenes, version='Old', 
                       species=TRUE, filter.no=1, rep.seq=NULL,
                       rff=FALSE, dep=NULL,
                       norm='TSS', level='OTU', intersect.no=4,
                       winsor=FALSE, winsor.qt=0.97,
                       ko.file=NULL, cog.file=NULL, ko.ann.file=NULL,
                       meta.sep='\t', quote="\"", comment="",
                       read.gg=FALSE,  seed=1234, ...) {
  # ko and cog file are not rarefied	
  # filter.no: filter the OTUs with read support less than filter.no (default is filtering singleton); singleton will not be filtered after rarefaction
  # winsorization and GMPR should be further studied. Current default is false and GMPR is on the genus level
  act.seq <- NULL
  set.seed(seed)
  
  meta.dat <- as.tibble(fread("maptest.txt"))
  
  # Load Tree
  tree.12 <- NULL
  if (!is.null(tree.file)) {
    tree.12 <- read.tree(tree.file)
    if (is.rooted(tree.12) == F) {
      tree.12 <- midpoint(tree.12)
    }
  }
  
  filter.no = 1
  biom <- read_biom(otu.file)
  otu.tab <- as(biom_data(biom), "matrix")
  otu.name <- observation_metadata(biom)
  
  otu.ind <- rowSums(otu.tab) > filter.no  # change otu.tab.12 != 0, rev:2016-06-20
  otu.tab <- otu.tab[otu.ind, ]
  # OTU names
  names(otu.name) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')
  otu.name.full <- as.matrix(otu.name[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')])
  otu.name.full <- otu.name.full[otu.ind, ]
  
  
  
  if (rff == TRUE) {
    cat('Rarefaction ...\n')
    if (is.null(dep)) {
      otu.tab.12 <- t(Rarefy(t(otu.tab.12))$otu.tab.rff)
    } else {
      otu.tab.12 <- t(Rarefy(t(otu.tab.12), dep)$otu.tab.rff)
    }	
    dep <- colSums(otu.tab.12)[1]
    cat("Depth ", dep, '\n')
    # Remove empty OTUs
    otu.ind <- rowSums(otu.tab.12) > 0  # rev:2016-06-28
    otu.tab.12 <- otu.tab.12[otu.ind, ]	
    otu.name.12 <- otu.name.12[otu.ind, ]
    otu.name.full <- otu.name.full[otu.ind, ]
    act.seq <- paste0(act.seq, 'R')
  } else {
    rff <- FALSE
    dep <- NULL
  }
  
  # Load mapping file
  if (load.map == TRUE) {
    samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))
    if (length(samIDs) == 0) {
      stop('Sample names in the meta file and biom are completely different?\n')
    } 
    if (length(samIDs) < length(colnames(otu.tab.12)) | length(samIDs) < length(rownames(meta.dat))) {
      warning('Sample names in the meta file and biom differ! May be due to rarefaction?\n')
    }
    meta.dat <- meta.dat[samIDs, ]
    otu.tab.12 <- otu.tab.12[, samIDs]
  } else {
    samIDs <- colnames(otu.tab.12)
  }
  
  # Create abundance list
  cat("Create taxa abundance list ...\n")
  abund.list.12 <- list()
  hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
  for (hierach in hierachs) {	
    if (hierach != 'Phylum') {
      single.names <- otu.name.12[, hierach]
      tax.family <- paste(otu.name.12[, 'Phylum'], single.names, sep=";")
      tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
    } else {
      tax.family <- otu.name.12[, 'Phylum']
      tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
    }
    family <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
    rownames(family) <- family[, 1]
    family <- as.matrix(family[, -1])
    abund.list.12[[hierach]] <- family
  }
  
  if (species) {
    abund.list.12[['Species']] <- otu.tab.12
    rownames(abund.list.12[['Species']]) <- paste0("OTU", rownames(otu.tab.12), ":", otu.name.12[, 'Phylum'], ";", otu.name.12[, 'Genus'])
  }
  
  cat('Normalize (size factor) ...\n')
  if (rff == TRUE) {
    cat('For rarefied data, the size factor for samples can still be different!\n')
  }
  
  if (level == 'OTU') {
    data <- otu.tab.12
  } else {
    if (level %in% names(abund.list.12)) {
      data <- abund.list.12[[level]]
    } else {
      data <- otu.tab.12
      level <-'OTU'
      warning('No or wrong level specified! OTU level will be used!\n')
    }
  }
  
  # Rarefaction/Normalizing factors are not calculated for functional data
  # Rev: 2017_01_19 the sample IDs for functional data are not ordered! Potential Danger! augment with NA
  if (!is.null(ko.file)) {
    
    cat("Load kegg file...\n")
    ko <- read_biom(ko.file)
    ko.dat <- as.matrix(biom_data(ko))
    
    if (sum(!(samIDs %in% colnames(ko.dat))) != 0) {
      missingIDs <- setdiff(samIDs, colnames(ko.dat))
      aug.mat <- matrix(NA, nrow(ko.dat), length(missingIDs))
      colnames(aug.mat) <- missingIDs
      rownames(aug.mat) <- rownames(ko.dat)
      ko.dat <- cbind(ko.dat, aug.mat) 
    }	
    
    ko.dat <- ko.dat[, samIDs]
    # Rarefaction?
    
    if (is.null(ko.ann.file)) {
      # Old - back compatability
      ko.ann <- observation_metadata(ko)
      ko.ann <- cbind(KEGG_Pathways1=sapply(ko.ann, function(x) x['KEGG_Pathways1']), 
                      KEGG_Pathways2=sapply(ko.ann, function(x) x['KEGG_Pathways2']), 
                      KEGG_Pathways3=sapply(ko.ann, function(x) x['KEGG_Pathways3']))
      rownames(ko.ann) <- rownames(ko.dat)		
      ko.ann[is.na(ko.ann)] <- 'Unclassified'
      
      hierachs <- c("KEGG_Pathways1", "KEGG_Pathways2", "KEGG_Pathways3")
      for (hierach in hierachs) {	
        tax.family <- ko.ann[, hierach]
        family <- aggregate(ko.dat, by=list(tax.family), FUN=sum)
        rownames(family) <- family[, 1]
        family <- as.matrix(family[, -1])
        abund.list.12[[hierach]] <- family
      }
    } else {
      # New
      load(ko.ann.file)
      #
      kos <- rownames(ko.dat)
      abund.list.12[["KEGG_Pathways3"]] <- NULL
      kos.id <- NULL
      for (ko.item in names(kegg.map)) {
        kos.common <- intersect(kos, kegg.map[[ko.item]])
        if (length(kos.common) != 0) {
          abund.list.12[["KEGG_Pathways3"]] <- rbind(abund.list.12[["KEGG_Pathways3"]], colSums(ko.dat[kos.common, , drop=F]))
          kos.id <- c(kos.id, ko.item)
        }
      }
      rownames(abund.list.12[["KEGG_Pathways3"]]) <- kos.id
      
      abund.list.12[["KEGG_Metabolism"]] <- abund.list.12[["KEGG_Pathways3"]][intersect(kos.id, unlist(kegg.ann[['Metabolism']])), ]
      rownames(abund.list.12[["KEGG_Metabolism"]]) <- paste0('M', rownames(abund.list.12[["KEGG_Metabolism"]])) 
      
      abund.list.12[["KEGG_Defense"]] <- NULL
      kos.id <- NULL
      for (ko.item in names(defense.map)) {
        kos.common <- intersect(kos, defense.map[[ko.item]])
        if (length(kos.common) != 0) {
          abund.list.12[["KEGG_Defense"]] <- rbind(abund.list.12[["KEGG_Defense"]], colSums(ko.dat[kos.common, , drop=F]))
          kos.id <- c(kos.id, ko.item)
        }
      }
      rownames(abund.list.12[["KEGG_Defense"]]) <- kos.id
      
      abund.list.12[["KEGG_Toxin"]] <- NULL
      kos.id <- NULL
      for (ko.item in names(toxin.map)) {
        kos.common <- intersect(kos, toxin.map[[ko.item]])
        if (length(kos.common) != 0) {
          abund.list.12[["KEGG_Toxin"]] <- rbind(abund.list.12[["KEGG_Toxin"]], colSums(ko.dat[kos.common, , drop=F]))
          kos.id <- c(kos.id, ko.item)
        }
      }
      rownames(abund.list.12[["KEGG_Toxin"]]) <- kos.id
    }
    
  }
  
  if (!is.null(cog.file)) {
    cat("Load cog file...\n")
    cog <- read_biom(cog.file)
    cog.dat <- as.matrix(biom_data(cog))
    
    if (sum(!(samIDs %in% colnames(cog.dat ))) != 0) {
      missingIDs <- setdiff(samIDs, colnames(cog.dat ))
      aug.mat <- matrix(NA, nrow(cog.dat), length(missingIDs))
      colnames(aug.mat) <- missingIDs
      rownames(aug.mat) <- rownames(cog.dat)
      cog.dat  <- cbind(cog.dat, aug.mat) 
    }	
    
    cog.dat <- cog.dat[, samIDs]
    
    # rarefaction?
    cog.ann <- observation_metadata(cog)
    hierachs <- c("COG_Category1", "COG_Category2")
    for (hierach in hierachs) {	
      tax.family <- sapply(cog.ann, function(x) x[hierach])
      family <- aggregate(cog.dat, by=list(tax.family), FUN=sum)
      rownames(family) <- family[, 1]
      family <- as.matrix(family[, -1])
      abund.list.12[[hierach]] <- family
    }
  }
  
  # Drop tree tips
  if (!is.null(tree.12)) {
    absent <- tree.12$tip.label[!(tree.12$tip.label %in% rownames(otu.tab.12))]
    if (length(absent) != 0) {
      tree.12 <- drop.tip(tree.12, absent)
      warning("The tree has OTUs not in the OTU table!")
    }
  }
  
  data.obj <- list(otu.tab=otu.tab.12, abund.list=abund.list.12, meta.dat=meta.dat, tree=tree.12,
                   otu.name=otu.name.12, otu.name.full=otu.name.full, 
                   size.factor=sf, norm.method=norm, norm.level=level, 
                   winsor=winsor, winsor.qt=winsor.qt,
                   rff=rff, rff.dep=dep, act.seq=act.seq,
                   call=match.call())
  
  data.obj <- list(otu.tab=otu.tab,otu.name=otu.name,otu.name.full=otu.name.full, tree=tree.12)
}

# Rev: 2016_09_22 Add load.map
# Rev: 2016_12_12 Reogranize rarefy, normalize, winsorize
load_data <- function (otu.file, map.file, tree.file=NULL,  load.map=TRUE, parseFunction=parse_taxonomy_greengenes, version='Old', 
                       species=TRUE, filter.no=1, rep.seq=NULL,
                       rff=FALSE, dep=NULL,
                       norm='TSS', level='OTU', intersect.no=4,
                       winsor=FALSE, winsor.qt=0.97,
                       ko.file=NULL, cog.file=NULL, ko.ann.file=NULL,
                       meta.sep='\t', quote="\"", comment="",
                       read.gg=FALSE,  seed=1234, ...) {
  # ko and cog file are not rarefied	
  # filter.no: filter the OTUs with read support less than filter.no (default is filtering singleton); singleton will not be filtered after rarefaction
  # winsorization and GMPR should be further studied. Current default is false and GMPR is on the genus level
  act.seq <- NULL
  set.seed(seed)
  if (load.map == TRUE) {
    cat("Load meta file...\n")
    if (grepl("csv$", map.file)) {
      meta.dat <- read.csv(map.file, header=T, check.names=F, row.names=1, comment=comment, quote=quote, ...)
    } else {
      meta.dat <- read.table(map.file, header=T, check.names=F, row.names=1, comment=comment, sep=meta.sep, quote=quote, ...)
    }
  } else {
    meta.dat <- NULL
  }
  
  # Load Tree
  if (!is.null(tree.file)) {
    cat("Load tree file ...\n")
    if (read.gg == F) {
      tree.12 <- read.tree(tree.file)
    } else {
      tree.12 <- read_tree_greengenes(tree.file)
    }
    
    if (is.rooted(tree.12) == F) {
      tree.12 <- midpoint(tree.12)
    }
    
  } else {
    tree.12 <- NULL
  }
  
  cat("Load OTU file...\n")  # Rewrite load new biom file, rev:2016-06-20
  if (version != 'New') {
    biom.obj <-  import_biom(otu.file, parseFunction = parseFunction)  		
    
    otu.tab.12 <- otu_table(biom.obj)@.Data
    otu.ind <- rowSums(otu.tab.12) > filter.no  # change otu.tab.12 != 0, rev:2016-06-20
    otu.tab.12 <- otu.tab.12[otu.ind, ]
    # OTU names
    otu.name.full <- as.matrix(biom.obj@tax_table[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')])
    otu.name.full <- otu.name.full[otu.ind, ]
    
    otu.name.12 <- otu.name.full[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
    otu.name.12[is.na(otu.name.12)] <- 'unclassified'
    otu.name.12 <- otu.name.12@.Data
    
    otu.name.12[, ] <- gsub('\\[', '', otu.name.12)
    otu.name.12[, ] <- gsub('\\]', '', otu.name.12)
    otu.name.full <- otu.name.full@.Data
  } else {
    temp <-  read_hdf5_biom(otu.file)
    otu.tab.12 <- matrix(unlist(temp$data), byrow=T, nrow=temp$shape[1], ncol=temp$shape[2])
    
    otu.ids <- sapply(temp$rows, function(x) x[['id']])
    sam.ids <- sapply(temp$columns, function(x) x[['id']])
    
    # From phyloseq: to be double checked!! Checked!
    if (all(sapply(sapply(temp$rows, function(i) {
      i$metadata
    }), is.null))) {
      otu.name.full <- NULL
    } else {
      taxlist = lapply(temp$rows, function(i) {
        parseFunction(i$metadata$taxonomy)
      })
      names(taxlist) = sapply(temp$rows, function(i) {
        i$id
      })
      otu.name.full = build_tax_table(taxlist)
    }
    
    otu.name.full <- otu.name.full@.Data
    
    rownames(otu.tab.12) <- otu.ids
    colnames(otu.tab.12) <- sam.ids
    
    otu.ind <- rowSums(otu.tab.12) > filter.no  # change otu.tab.12 != 0, rev:2016-06-20
    otu.tab.12 <- otu.tab.12[otu.ind, ]	
    otu.name.full <- otu.name.full[otu.ind, ]
    otu.name.12 <- otu.name.full[, c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')]
    otu.name.12[is.na(otu.name.12)] <- 'unclassified'
    otu.name.12[, ] <- gsub('\\[', '', otu.name.12)
    otu.name.12[, ] <- gsub('\\]', '', otu.name.12)
  }
  
  if (rff == TRUE) {
    cat('Rarefaction ...\n')
    if (is.null(dep)) {
      otu.tab.12 <- t(Rarefy(t(otu.tab.12))$otu.tab.rff)
    } else {
      otu.tab.12 <- t(Rarefy(t(otu.tab.12), dep)$otu.tab.rff)
    }	
    dep <- colSums(otu.tab.12)[1]
    cat("Depth ", dep, '\n')
    # Remove empty OTUs
    otu.ind <- rowSums(otu.tab.12) > 0  # rev:2016-06-28
    otu.tab.12 <- otu.tab.12[otu.ind, ]	
    otu.name.12 <- otu.name.12[otu.ind, ]
    otu.name.full <- otu.name.full[otu.ind, ]
    act.seq <- paste0(act.seq, 'R')
  } else {
    rff <- FALSE
    dep <- NULL
  }
  
  # Load mapping file
  if (load.map == TRUE) {
    samIDs <- intersect(rownames(meta.dat), colnames(otu.tab.12))
    if (length(samIDs) == 0) {
      stop('Sample names in the meta file and biom are completely different?\n')
    } 
    if (length(samIDs) < length(colnames(otu.tab.12)) | length(samIDs) < length(rownames(meta.dat))) {
      warning('Sample names in the meta file and biom differ! May be due to rarefaction?\n')
    }
    meta.dat <- meta.dat[samIDs, ]
    otu.tab.12 <- otu.tab.12[, samIDs]
  } else {
    samIDs <- colnames(otu.tab.12)
  }
  
  # Create abundance list
  cat("Create taxa abundance list ...\n")
  abund.list.12 <- list()
  hierachs <- c('Phylum', 'Class', 'Order', 'Family', 'Genus')
  for (hierach in hierachs) {	
    if (hierach != 'Phylum') {
      single.names <- otu.name.12[, hierach]
      tax.family <- paste(otu.name.12[, 'Phylum'], single.names, sep=";")
      tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- paste0('Unclassified_', hierach)
    } else {
      tax.family <- otu.name.12[, 'Phylum']
      tax.family[grepl('unclassified', tax.family, ignore.case=T)] <- 'Unclassified_Phylum'
    }
    family <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
    rownames(family) <- family[, 1]
    family <- as.matrix(family[, -1])
    abund.list.12[[hierach]] <- family
  }
  
  if (species) {
    abund.list.12[['Species']] <- otu.tab.12
    rownames(abund.list.12[['Species']]) <- paste0("OTU", rownames(otu.tab.12), ":", otu.name.12[, 'Phylum'], ";", otu.name.12[, 'Genus'])
  }
  
  cat('Normalize (size factor) ...\n')
  if (rff == TRUE) {
    cat('For rarefied data, the size factor for samples can still be different!\n')
  }
  
  if (level == 'OTU') {
    data <- otu.tab.12
  } else {
    if (level %in% names(abund.list.12)) {
      data <- abund.list.12[[level]]
    } else {
      data <- otu.tab.12
      level <-'OTU'
      warning('No or wrong level specified! OTU level will be used!\n')
    }
  }
  
  if (norm == 'GMPR') {
    sf <- GMPR(data, intersect.no=intersect.no)
    warning('GMPR is only suitable for samples from the same body location!\n')
    norm <- 'GMPR'
    names(sf) <- colnames(data)
    act.seq <- paste0(act.seq, 'N')
  } else {
    if (norm == 'TSS') {
      sf <- colSums(data)
      norm <- 'TSS'
      act.seq <- paste0(act.seq, 'N')
    } else {
      warning('Normalization method not specified or unknown! TSS is used!\n')
      sf <- colSums(data)
      norm <- 'TSS'
      act.seq <- paste0(act.seq, 'N')
    }
    
  }
  
  if (winsor == TRUE) {
    act.seq <- paste0(act.seq, 'W')
    if (rff == TRUE) {
      warning('Winsorization after rarefaction will make the data have different total numbers!\n')
    }
    cat('Winsorize ...\n')
    if (is.null(winsor.qt)) {
      winsor.qt <- 0.97
    }
    # Addressing the outlier (97% percent) or at least one outlier
    abund.list.12 <- sapply(abund.list.12, function(genus) {
      genus.p <- t(t(genus) / sf)
      genus.p <- apply(genus.p, 1, function(x) {
        cutoff <- quantile(x, winsor.qt)
        x[x >= cutoff] <- cutoff
        x
      }
      )
      # column/row switch
      genus.w <- t(round(genus.p * sf))
    })
    # OTU table
    otu.tab.12.p <- t(t(otu.tab.12) / sf)
    otu.tab.12.p <- apply(otu.tab.12.p, 1, function(x) {
      cutoff <- quantile(x, winsor.qt)
      x[x >= cutoff] <- cutoff
      x
    }
    )
    # column/row switch
    otu.tab.12 <- t(round(otu.tab.12.p * sf))
  } else {
    winsor <- FALSE
    winsor.qt <- NULL
  }
  
  # Rarefaction/Normalizing factors are not calculated for functional data
  # Rev: 2017_01_19 the sample IDs for functional data are not ordered! Potential Danger! augment with NA
  if (!is.null(ko.file)) {
    cat("Load kegg file...\n")
    ko <- read_biom(ko.file)
    ko.dat <- as.matrix(biom_data(ko))
    
    if (sum(!(samIDs %in% colnames(ko.dat))) != 0) {
      missingIDs <- setdiff(samIDs, colnames(ko.dat))
      aug.mat <- matrix(NA, nrow(ko.dat), length(missingIDs))
      colnames(aug.mat) <- missingIDs
      rownames(aug.mat) <- rownames(ko.dat)
      ko.dat <- cbind(ko.dat, aug.mat) 
    }
    
    ko.dat <- ko.dat[, samIDs]
    # Rarefaction?
    
    if (is.null(ko.ann.file)) {
      # Old - back compatability
      ko.ann <- observation_metadata(ko)
      ko.ann <- cbind(KEGG_Pathways1=sapply(ko.ann, function(x) x['KEGG_Pathways1']), 
                      KEGG_Pathways2=sapply(ko.ann, function(x) x['KEGG_Pathways2']), 
                      KEGG_Pathways3=sapply(ko.ann, function(x) x['KEGG_Pathways3']))
      rownames(ko.ann) <- rownames(ko.dat)		
      ko.ann[is.na(ko.ann)] <- 'Unclassified'
      
      hierachs <- c("KEGG_Pathways1", "KEGG_Pathways2", "KEGG_Pathways3")
      for (hierach in hierachs) {	
        tax.family <- ko.ann[, hierach]
        family <- aggregate(ko.dat, by=list(tax.family), FUN=sum)
        rownames(family) <- family[, 1]
        family <- as.matrix(family[, -1])
        abund.list.12[[hierach]] <- family
      }
    } else {
      # New
      load(ko.ann.file)
      #
      kos <- rownames(ko.dat)
      abund.list.12[["KEGG_Pathways3"]] <- NULL
      kos.id <- NULL
      for (ko.item in names(kegg.map)) {
        kos.common <- intersect(kos, kegg.map[[ko.item]])
        if (length(kos.common) != 0) {
          abund.list.12[["KEGG_Pathways3"]] <- rbind(abund.list.12[["KEGG_Pathways3"]], colSums(ko.dat[kos.common, , drop=F]))
          kos.id <- c(kos.id, ko.item)
        }
      }
      rownames(abund.list.12[["KEGG_Pathways3"]]) <- kos.id
      
      abund.list.12[["KEGG_Metabolism"]] <- abund.list.12[["KEGG_Pathways3"]][intersect(kos.id, unlist(kegg.ann[['Metabolism']])), ]
      rownames(abund.list.12[["KEGG_Metabolism"]]) <- paste0('M', rownames(abund.list.12[["KEGG_Metabolism"]])) 
      
      abund.list.12[["KEGG_Defense"]] <- NULL
      kos.id <- NULL
      for (ko.item in names(defense.map)) {
        kos.common <- intersect(kos, defense.map[[ko.item]])
        if (length(kos.common) != 0) {
          abund.list.12[["KEGG_Defense"]] <- rbind(abund.list.12[["KEGG_Defense"]], colSums(ko.dat[kos.common, , drop=F]))
          kos.id <- c(kos.id, ko.item)
        }
      }
      rownames(abund.list.12[["KEGG_Defense"]]) <- kos.id
      
      abund.list.12[["KEGG_Toxin"]] <- NULL
      kos.id <- NULL
      for (ko.item in names(toxin.map)) {
        kos.common <- intersect(kos, toxin.map[[ko.item]])
        if (length(kos.common) != 0) {
          abund.list.12[["KEGG_Toxin"]] <- rbind(abund.list.12[["KEGG_Toxin"]], colSums(ko.dat[kos.common, , drop=F]))
          kos.id <- c(kos.id, ko.item)
        }
      }
      rownames(abund.list.12[["KEGG_Toxin"]]) <- kos.id
    }
    
  }
  
  if (!is.null(cog.file)) {
    cat("Load cog file...\n")
    cog <- read_biom(cog.file)
    cog.dat <- as.matrix(biom_data(cog))
    if (sum(!(samIDs %in% colnames(cog.dat ))) != 0) {
      missingIDs <- setdiff(samIDs, colnames(cog.dat ))
      aug.mat <- matrix(NA, nrow(cog.dat), length(missingIDs))
      colnames(aug.mat) <- missingIDs
      rownames(aug.mat) <- rownames(cog.dat)
      cog.dat  <- cbind(cog.dat, aug.mat) 
    }	
    
    cog.dat <- cog.dat[, samIDs]
    # rarefaction?
    cog.ann <- observation_metadata(cog)
    hierachs <- c("COG_Category1", "COG_Category2")
    for (hierach in hierachs) {	
      tax.family <- sapply(cog.ann, function(x) x[hierach])
      family <- aggregate(cog.dat, by=list(tax.family), FUN=sum)
      rownames(family) <- family[, 1]
      family <- as.matrix(family[, -1])
      abund.list.12[[hierach]] <- family
    }
  }
  
  # Drop tree tips
  if (!is.null(tree.12)) {
    absent <- tree.12$tip.label[!(tree.12$tip.label %in% rownames(otu.tab.12))]
    if (length(absent) != 0) {
      tree.12 <- drop.tip(tree.12, absent)
      warning("The tree has OTUs not in the OTU table!")
    }
  }
  
  data.obj <- list(otu.tab=otu.tab.12, abund.list=abund.list.12, meta.dat=meta.dat, tree=tree.12,
                   otu.name=otu.name.12, otu.name.full=otu.name.full, 
                   size.factor=sf, norm.method=norm, norm.level=level, 
                   winsor=winsor, winsor.qt=winsor.qt,
                   rff=rff, rff.dep=dep, act.seq=act.seq,
                   call=match.call())
}

# Rev: 2016_11_28, Bray-curtis use normalized data.
construct_distance <- function (data.obj, unifrac.file=NULL,  Phylum='All', dist.RData=NULL, save.RData=NULL, 
                                filter.no=0, rff=FALSE, dep=NULL, seed=1234) {
  set.seed(seed)
  
  if (!is.null(dist.RData)) {
    load(dist.RData, envir=.GlobalEnv )
  } else {
    dist.list.12 <- list()
    cat("Generalized UniFrac ...\n")
    
    otu.tab <- t(data.obj$otu.tab)
    
    if (rff == TRUE) {
      if (is.null(dep)) {
        otu.tab <- Rarefy(otu.tab)$otu.tab.rff
      } else {
        otu.tab <- Rarefy(otu.tab, dep)$otu.tab.rff
      }	
    }
    
    if (Phylum != 'All') {
      ind <- data.obj$otu.name[, 'Phylum'] == Phylum
      otu.tab <- otu.tab[, ind]
    }
    
    # Filter otus with reads <= filter.no
    otu.tab <- otu.tab[, colSums(otu.tab) > filter.no]
    
    # Remove samples with no reads
    if (sum(rowSums(otu.tab) == 0) >= 1) {
      otu.tab <- otu.tab[rowSums(otu.tab) != 0, ]
      warning('Some samples do not have reads after rarefaction! Please be careful!\n')
    }
    
    # To make sure the OTUs in otu.tab are in the tree (rev:2016-06-19)
    
    if (sum(!(colnames(otu.tab) %in% data.obj$tree$tip.label))) {
      warning('Some OTU names are not in the tree! An intersection set will be used!\n')	
    }
    common.otus <- intersect(colnames(otu.tab), data.obj$tree$tip.label)
    unifrac12 <- GUniFrac(otu.tab[, common.otus], data.obj$tree)$unifracs
    
    dist.list.12[['WUniFrac']] <- unifrac12[, , 'd_1']
    dist.list.12[['GUniFrac']] <- unifrac12[, , 'd_0.5']
    if (is.null(unifrac.file)) {
      dist.list.12[['UniFrac']] <- unifrac12[, , 'd_UW']
    } else {
      # The orders may be different
      dist.list.12[['UniFrac']] <- as.matrix(read.table(unifrac.file, row.names=1, header=T)) # Rarefaction
    }
    
    # Need speed up
    # Suggest using rarefied counts 
    # If case/control has different sequencing depth, it will result in false clustering! 
    # Rev: 2016_11_28 Use normalized data for BC distance calculation to reduce noise in case of variable library size
    dist.list.12[['BC']] <-as.matrix(vegdist(otu.tab / rowSums(otu.tab)))
    
    genus <- t(data.obj$abund.list[['Genus']])
    genus <- genus / rowSums(genus)
    dist.list.12[['Euc']] <- as.matrix(dist(genus))
    
    genus <- sqrt(genus)
    dist.list.12[['Hel']] <-as.matrix(dist(genus))
    #		dist.list.12[['JS']] <- as.matrix(distance(otu_table(data.obj$abund.list[['Genus']], taxa_are_rows=T), method='jsd'))
    
    if (!is.null(save.RData)) {
      save(dist.list.12, file=save.RData)
    }
  }
  
  return(dist.list.12)
}


is.na.null <- function (x) {
  if (is.null(x)) {
    return(TRUE)
  } else {
    if (is.na(x)[1]) {
      return(TRUE)
    }  else {
      return(FALSE)
    }
  }
  
}

# Rev: 2016_09_26 remove empty OTUs/taxa
# Rev: 2016_12_01 add more logical controls
subset_data <- function (data.obj, samIDs) {
  
  # Rev: 2016_1_19 to add error protection
  # Transform logical samIDs into characer samIDs
  if (is.logical(samIDs) | is.numeric(samIDs)) {
    samIDs <- rownames(data.obj$meta.dat)[samIDs]
  }
  
  data.obj$meta.dat <- data.obj$meta.dat[samIDs, , drop=FALSE]
  
  if (!is.na.null(data.obj$otu.tab)) {
    
    data.obj$otu.tab <- data.obj$otu.tab[, samIDs, drop=FALSE]
    data.obj$otu.tab <- data.obj$otu.tab[rowSums(data.obj$otu.tab) != 0, , drop=FALSE]
    data.obj$otu.name <- data.obj$otu.name[rownames(data.obj$otu.tab), , drop=FALSE]
    
    if (!is.na.null(data.obj$otu.name.full)) {
      data.obj$otu.name.full <- data.obj$otu.name.full[rownames(data.obj$otu.tab), , drop=FALSE]
    }
  }
  
  if (!is.na.null(data.obj$abund.list)) {
    data.obj$abund.list <- lapply(data.obj$abund.list, function(x) {
      xx <- x[, samIDs, drop=FALSE]
      xx <- xx[rowSums(xx) != 0, , drop=FALSE]
    })
  }
  
  if (!is.na.null(data.obj$size.factor)) {
    data.obj$size.factor <- data.obj$size.factor[samIDs]
  }
  
  if (!is.na.null(data.obj$ko.list)) {
    data.obj$ko.list <- lapply(data.obj$ko.list, function(x) {
      xx <- x[, samIDs, drop=FALSE]
      xx <- xx[rowSums(xx) != 0, , drop=FALSE]
    })
  }
  
  if (!is.na.null(data.obj$cog.list)) {
    data.obj$cog.list <- lapply(data.obj$cog.list, function(x) {
      xx <- x[, samIDs, drop=FALSE]
      xx <- xx[rowSums(xx) != 0, , drop=FALSE]
    })
  }
  return(data.obj)
}

# Rev: 2016_12_01 add more logical controls
# Rev: 2016_02_01 fix one error
subset_dist <- function (dist.obj, samIDs) {
  
  # Rev: 2016_1_19 to add error protection
  # Transform logical samIDs into character samIDs
  if (is.logical(samIDs) | is.numeric(samIDs)) {
    samIDs <- rownames(dist.obj[[1]])[samIDs]
  }
  
  lapply(dist.obj, function(x) {
    if(!is.na.null(x)){
      x <- x[samIDs, samIDs]
    } else {
      x
    }
    x
  })
}


perform_sequence_stat_analysis <- function (data.obj, ann='') {
  sink(paste0('Sequence_Analysis_Statistics_', ann, '.txt'))
  otu.tab <- data.obj$otu.tab
  
  # Sequencing depth
  otu.abund <- rowSums(otu.tab)
  sam.abund <- colSums(otu.tab)
  otu.prev <- rowSums(otu.tab!=0)/ncol(otu.tab)
  
  otu.abund <- otu.abund[otu.abund >= 1]
  sam.abund <- sam.abund[sam.abund >= 1]
  cat('This data set contains ', length(sam.abund), ' samples after quality controls.\n')
  cat('16S rDNA targeted sequencing yields ', mean(sam.abund), 'reads/sample on average (range:', min(sam.abund), '-', max(sam.abund), ').')
  cat('Clustering of these 16S sequence tags produces ', sum(otu.abund > 0), ' OTUs at 97% similarity level.')
  
  
  
  phy.abund <- data.obj$abund.list[['Phylum']]
  fam.abund <- data.obj$abund.list[['Family']]
  gen.abund <- data.obj$abund.list[['Genus']]
  
  phy.prev <- rowSums(phy.abund != 0) / ncol(phy.abund)
  fam.prev <- rowSums(fam.abund != 0) / ncol(phy.abund)
  gen.prev <- rowSums(gen.abund != 0) / ncol(phy.abund)
  
  phy.abund <- rowMeans(t(t(phy.abund) / sam.abund))
  fam.abund <- rowMeans(t(t(fam.abund) / sam.abund))
  gen.abund <- rowMeans(t(t(gen.abund) / sam.abund))
  
  cat('These OTUs belong to ', sum(phy.abund > 0), ' phyla,', sum(fam.abund > 0), ' families and ', sum(gen.abund > 0), 'genera.\n\n')
  
  phy.prev <- sort(phy.prev, decr=T)
  phy.prev <- round(phy.prev[phy.prev >= 0.05] * 100, 2)
  
  fam.prev <- sort(fam.prev, decr=T)
  fam.prev <- round(fam.prev[fam.prev >= 0.05] * 100, 2)
  
  gen.prev <- sort(gen.prev, decr=T)
  gen.prev <- round(gen.prev[gen.prev >= 0.05] * 100, 2)
  
  # Rev: 2017_02_19 ' ' -> '\n'
  cat('\nThe most prevalent phyla are:\n', paste(paste0(names(phy.prev), '(', phy.prev, '%)'), collapse='\n'), '\n')
  cat('\nThe most prevalent families are:\n', paste(paste0(names(fam.prev), '(', fam.prev, '%)'), collapse='\n'), '\n')
  cat('\nand the most prevalent genera are:\n', paste(paste0(names(gen.prev), '(', gen.prev, '%)'), collapse='\n'), '\n\n')
  
  phy.abund <- sort(phy.abund, decr=T)
  phy.abund <- round(phy.abund[phy.abund >= 0.05] * 100, 2)
  
  fam.abund <- sort(fam.abund, decr=T)
  fam.abund <- round(fam.abund[fam.abund >= 0.05] * 100, 2)
  
  gen.abund <- sort(gen.abund, decr=T)
  gen.abund <- round(gen.abund[gen.abund >= 0.05] * 100, 2)
  
  cat('\nThe most abundant phyla are ', paste(paste0(names(phy.abund), '(', phy.abund, '%)'), collapse=' '), ';')
  cat('\nThe most abundant families are ', paste(paste0(names(fam.abund), '(', fam.abund, '%)'), collapse=' '), ';')
  cat('\nand the most abundant genera are ', paste(paste0(names(gen.abund), '(', gen.abund, '%)'), collapse=' '), '.')
  sink()
  sink(paste0('Sequence_Analysis_Statistics_table_', ann, '.tsv'))
  write.table(cbind(read.table(text=names(phy.prev)), unname(phy.prev), "Phylum", "Prevalence"), row.names=FALSE, col.names=FALSE)
  write.table(cbind(read.table(text=names(fam.prev)), unname(fam.prev), "Family", "Prevalence"), row.names=FALSE, col.names=FALSE)
  write.table(cbind(read.table(text=names(gen.prev)), unname(gen.prev), "Genus", "Prevalence"), row.names=FALSE, col.names=FALSE)
  write.table(cbind(read.table(text=names(phy.abund)), unname(phy.abund), "Phylum", "Abundance"), row.names=FALSE, col.names=FALSE)
  write.table(cbind(read.table(text=names(fam.abund)), unname(fam.abund), "Family", "Abundance"), row.names=FALSE, col.names=FALSE)
  write.table(cbind(read.table(text=names(gen.abund)), unname(gen.abund), "Genus", "Abundance"), row.names=FALSE, col.names=FALSE)
  sink()
  
  png(paste0('Sequence_Analysis_Statistics_', ann, '.png'), height=600, width=900)
  
  obj1 <- ggplot2::ggplot(data=data.frame(x=sam.abund), aes(x=x)) + geom_histogram(col='black', fill='gray')  + ylab('Frequency') + xlab('Sequencing depth') + theme_bw()
  
  otu.tab <- data.obj$otu.tab
  map <- data.obj$meta.dat
  colnames(otu.tab) <-  map[[ann]]
  df <- data.frame(Group=names(colSums(otu.tab)), coverage=colSums(otu.tab))
  
  obj2 <- ggplot2::ggplot(df, aes(x=Group, y=log10(coverage), col=Group)) + geom_boxplot(position=position_dodge(width=0.75), outlier.colour = NA) +
    geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) + theme_bw()
  multiplot(obj1,obj2, cols=1)
  dev.off()
  sink(paste0('Sequence_Analysis_Statistics_log10Coverage_Association_', ann, '.tsv'))
  map$log10_coverage <- log10(colSums(data.obj$otu.tab))
  lm.obj <- lm(as.formula(paste('log10_coverage ~ ', ann)), map)
  prmatrix(summary(lm.obj)$coefficients)
  sink()
}

# Rev: 2017_08_23 add automatically create phylo.obj
generate_rarefy_curve <- function (data.obj, phylo.obj=NULL, grp.name, depth=NULL, npoint=10, iter.no=5,
                                   measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), ann='', gg.cmd="theme(legend.justification=c(1,0), legend.position=c(1,0))", wid=5, hei=5) {
  cat("Create rarefaction curves!\n")
  # Rev: 2017_08_23
  if (is.null(phylo.obj)) {
    phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
                          tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
  }
  if (is.null(depth)) {
    depth <- min(sample_sums(phylo.obj))
    phylo.even <- rarefy_even_depth(phylo.obj, rngseed=12345)
  } else {
    if (depth > min(sample_sums(phylo.obj))) {
      ind <- sample_sums(phylo.obj) >= depth
      cat(sum(!ind), " samples do not have sufficient number of reads!\n")
      sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ] 
      data.obj <- subset_data(data.obj, ind)
    }
    phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
  }
  
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  
  # Rev: 2016_12_12
  if (is.character(grp)) {
    grp <- factor(grp)
  }
  
  if (!is.factor(grp)) {
    stop('Rarefaction curve needs a factor!\n')
  }
  res <- NULL
  incr <- depth %/% npoint
  sink('temp.txt')
  for (dep in c(10, incr*(1:npoint))) {
    x <- 0
    for (i in 1:iter.no) {
      phylo.even <- rarefy_even_depth(phylo.obj, dep, rngseed=12345+i)
      x <- x + estimate_richness(phylo.even, measures=measures)
    }
    
    res <- rbind(res, t(x[, measures, drop=F]/iter.no))
  }
  colnames(res) <- rownames(df)
  sink()
  res_list <- list()
  
  for (i in 1:length(measures)) {
    measure <- measures[i]
    cat("Measure: ", measure, "\n")
    res2 <- res[(0:(npoint))*length(measures)+i, , drop=F]
    m <- t(apply(res2, 1, function(x) tapply(x, grp, mean)))
    se <- t(apply(res2, 1, function(x) tapply(x, grp, function(y) sd(y)/sqrt(length(y)))))
    uci <- m+se
    lci <- m-se
    
    m <- melt(m)
    uci <- melt(uci)
    lci <- melt(lci)
    
    res2 <- cbind(c(10, incr*(1:npoint)), m[, 2:3], uci[, 3], lci[, 3])
    colnames(res2) <- c('Depth', 'Group', 'mean', 'max', 'min')
    
    res2 <- as.data.frame(res2)
    res2$Group <- factor(res2$Group, levels=levels(grp))

    res_list[[measure]] <- res2
    res2
    
  }
  mres_list = melt(res_list, measure.vars=c("mean", "max", "min"))
  cres_list <- dcast(mres_list, L1 + Group + Depth ~ variable, fun.aggregate=mean)
  obj <- ggplot(cres_list, aes(x=Depth, y=mean, color=Group, group=Group)) +
    geom_errorbar(aes(ymin=min, ymax=max), alpha=0.5, width=.25, position=position_dodge(.2)) + 
    geom_line() + 
    geom_point(size=3, shape=21, fill="white") +
    labs(y="Alpha diversity") +
    facet_wrap(~ L1, scale="free_y") + theme_bw()
  return(obj)
}

# Rev: 2016_09_10
# Rev: 2016_11_28
# Rev: 2017_04_18
generate_alpha_boxplot <- function (data.obj, phylo.obj=NULL, rarefy=TRUE, depth=NULL, grp.name, strata=NULL, 
                                    measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), gg.cmd=NULL, ann='', subject=NULL, p.size=2.5, l.size=0.5,
                                    hei = NULL, wid = NULL) {	
  # Rev: 2017_08_23
  if (is.null(phylo.obj)) {
    phylo.obj <- phyloseq(otu_table(data.obj$otu.tab, taxa_are_rows=T), phy_tree(data.obj$tree), 
                          tax_table(data.obj$otu.name), sample_data(data.obj$meta.dat))
  }
  # To be completed - jetter when strata is not null
  if (rarefy == TRUE) {
    if (is.null(depth)) {
      depth <- min(sample_sums(phylo.obj))
    } else {
      if (depth > min(sample_sums(phylo.obj))) {
        ind <- (sample_sums(phylo.obj) >= depth)
        cat(sum(!ind), " samples do not have sufficient number of reads!\n")
        
        sample_data(phylo.obj) <- sample_data(phylo.obj)[ind, ]
        data.obj <- subset_data(data.obj, ind)
      }
    }
    
    phylo.even <- rarefy_even_depth(phylo.obj, depth, rngseed=12345)
    est_rich <- estimate_richness(phylo.even, measures=measures)
  } else {
    est_rich <- estimate_richness(phylo.obj, measures=measures)
  } 
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  obj_list <- list()
  df = data.frame(Value=est_rich[, measures], Group=grp)
  mdf <- melt(df, id="Group")
  obj <- ggplot(mdf, aes(x=as.factor(Group), y=value, col=as.factor(Group), group=as.factor(Group))) + 
    geom_boxplot(position=position_dodge(width=0.75), outlier.colour = NA) + 
    geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) + 
    labs(y="Alpha Diversity", x="Group") + 
    facet_wrap(~ variable, scale="free_y") + 
    theme_bw()
  return(obj)
}

# New: 2018_03_09 Removing phyloseq dependency and add 'subject' parameter, and remove 'model' parameter 
perform_alpha_test2 <- function (data.obj, alpha.obj=NULL, rarefy=TRUE, depth=NULL, iter.no=5, 
                                 measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'),  model='lm', 
                                 formula=NULL, grp.name=NULL, adj.name=NULL, subject=NULL, ann='', seed=123, ...) {
  if (is.null(alpha.obj)) {
    alpha.obj <- generate_alpha_diversity(data.obj,  rarefy=rarefy, depth=depth, iter.no=iter.no, measures=measures, seed=seed)
    
  } else {
    if (sum(!(rownames(alpha.obj) %in% rownames(data.obj$meta.dat))) != 0){
      stop("alpha.obj contains samples not in data.obj!\n")
    } 
  }
  
  
  x <- alpha.obj
  df <- data.obj$meta.dat[rownames(alpha.obj), ]
  
  if (is.null(subject)) {
    model <- 'lm'
  } else {
    model <- 'lme'
  }
  
  result <- list()
  fitted.obj <- list()
  
  
  # variables to adjust always come first in anova analyses
  if (is.null(formula)) {
    if (is.null(adj.name)) {
      formula <- paste('~', grp.name)
    } else {
      formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
    }
  }
  
  for (measure in measures) {
    cat("Alpha diversity:", measure, "\n")
    xx <- x[, measure]
    if (model == 'lm') {
      cat('Linear model:\n')
      
      lm.obj <- lm(as.formula(paste('xx ', formula)), df, ...)
      prmatrix(summary(lm.obj)$coefficients)
      cat('\nANOVA:\n')
      print(anova(lm.obj))
      cat('\n')
      fitted.obj[[measure]] <- lm.obj
    }
    if (model == 'lme') {
      df$xx <- xx
      cat('Linear mixed effects model:\n')
      lm.obj <- lme(as.formula(paste('xx ', formula)), df, method='ML', random=as.formula(paste0(' ~ 1 | ', subject)), ...)
      prmatrix(summary(lm.obj)$tTable)
      cat('\nANOVA:\n')
      print(anova(lm.obj))
      cat('\n')
      fitted.obj[[measure]] <- lm.obj
    }
    cat("\n")
  }
  
  result$fitted.obj <- fitted.obj
  result$alpha.diversity <- x
  # Rev: 2016_12_25
  return(result)
  
}

# New: 2018_02_25
estimate_richness_ <- function (OTU,  measures = NULL) {
  # Modified over Phyloseq
  
  renamevec = c("Observed", "Chao1", "ACE", "Shannon", "Simpson", 
                "InvSimpson", "Fisher")
  names(renamevec) <- c("S.obs", "S.chao1", "S.ACE", "shannon", 
                        "simpson", "invsimpson", "fisher")
  if (is.null(measures)) {
    measures = as.character(renamevec)
  }
  if (any(measures %in% names(renamevec))) {
    measures[measures %in% names(renamevec)] <- renamevec[names(renamevec) %in% 
                                                            measures]
  }
  if (!any(measures %in% renamevec)) {
    stop("None of the `measures` you provided are supported. Try default `NULL` instead.")
  }
  outlist = vector("list")
  estimRmeas = c("Chao1", "Observed", "ACE")
  if (any(estimRmeas %in% measures)) {
    outlist <- c(outlist, list(t(data.frame(vegan::estimateR(OTU)))))
  }
  if ("Shannon" %in% measures) {
    outlist <- c(outlist, list(shannon = vegan::diversity(OTU, index = "shannon")))
  }
  if ("Simpson" %in% measures) {
    outlist <- c(outlist, list(simpson = vegan::diversity(OTU, index = "simpson")))
  }
  if ("InvSimpson" %in% measures) {
    outlist <- c(outlist, list(invsimpson = vegan::diversity(OTU, 
                                                             index = "invsimpson")))
  }
  if ("Fisher" %in% measures) {
    fisher = tryCatch(vegan::fisher.alpha(OTU, se = TRUE), warning = function(w) {
      warning("estimate_richness: Warning in fisher.alpha(). See `?fisher.fit` or ?`fisher.alpha`. Treat fisher results with caution")
      suppressWarnings(vegan::fisher.alpha(OTU, se = TRUE)[, c("alpha", "se")])
    })
    if (!is.null(dim(fisher))) {
      colnames(fisher)[1:2] <- c("Fisher", "se.fisher")
      outlist <- c(outlist, list(fisher))
    }
    else {
      outlist <- c(outlist, Fisher = list(fisher))
    }
  }
  out = do.call("cbind", outlist)
  namechange = intersect(colnames(out), names(renamevec))
  colnames(out)[colnames(out) %in% namechange] <- renamevec[namechange]
  colkeep = sapply(paste0("(se\\.){0,}", measures), grep, colnames(out), 
                   ignore.case = TRUE)
  out = out[, sort(unique(unlist(colkeep))), drop = FALSE]
  out <- as.data.frame(out)
  return(out)
}

# New: 2018_02_25
generate_alpha_diversity <- function (data.obj,  rarefy=TRUE, depth=NULL, iter.no=5,
                                      measures=c('Observed', 'Chao1', 'Shannon', 'InvSimpson'), seed=123) {	
  # Rev: 2017_08_23
  
  OTU <- t(data.obj$otu.tab)
  depths <- rowSums(OTU)
  result <- list()
  if (rarefy == TRUE) {
    if (is.null(depth)) {
      depth <- min(depths)
    } else {
      if (depth > min(depths)) {
        ind <- depths >= depth
        cat(sum(!ind), " samples do not have sufficient number of reads!\n")
        # data.obj <- subset_data(data.obj, ind)
      }
    }
    
    set.seed(123)
    x <- 0 
    for (i in 1:iter.no) {
      OTU2 <- Rarefy(OTU, depth = depth)$otu.tab.rff
      x <- x + estimate_richness_(OTU2, measures=measures)
    }
    x <- x / iter.no
    rownames(x) <- rownames(OTU2)
  } else {
    x <- estimate_richness_(OTU, measures=measures)
    rownames(x) <- rownames(OTU)
  }
  
  return(x)
}

#Note: this function is NOT general purpose like generate_ordination()
generate_ordination2 <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
                                  grp.name, adj.name=NULL, emp.lev=NULL, strata=NULL, pc.separate, pca.method='cmd', ann=NULL, sub=NULL,
                                  clab=1.0, cex.pt=1.25, ellipse=T, cstar= 1, wid=900, hei=600, ...) {
  # Implment strata
  # To be completed, add continuous case
  strata0 <- strata
  
  df <- data.obj$meta.dat
  grp <- factor(df[, grp.name])
  
  if (is.null(emp.lev)) {
    grp <- factor(grp, levels=c(setdiff(levels(grp), emp.lev), emp.lev))
  }
  
  if (!is.null(strata)) {
    strata <- factor(df[, strata])
  } else {
    strata <- factor(grp)
  }
  y <- list()
  for (dist.name in dist.names) {
    dist.temp <- dist.obj[[dist.name]]
    if (!is.null(adj.name)) {
      adj <- as.data.frame(df[, adj.name])
      obj <- cmdscale(as.dist(dist.temp), k=ncol(dist.temp)-1)
      dat2 <- apply(obj, 2, function(x) resid(lm(x ~ ., data=adj)))
      dist.temp <- dist(dat2)
    } 
    if (pca.method == 'cmd') {
      obj <- cmdscale(as.dist(dist.temp), k=2, eig=T)
      pve <- round(obj$eig[1:2]/sum(abs(obj$eig))*100, 1)
      y[[dist.name]] <- cbind(obj$points[, 1], obj$points[, 2])
      xlab <- paste0('PC1(', pve[1], '%)')
      ylab <- paste0('PC2(', pve[2], '%)')
    } 
    
    if (pca.method == 'nmds') {
      obj <- metaMDS(as.dist(dist.temp), k=2)
      y[[dist.name]] <- cbind(obj$points[, 1], obj$points[, 2])
      xlab <- 'NMDS1'
      ylab <- 'NMDS2'
    } 
    
    colnames(y[[dist.name]]) <- c("PC1", "PC2")
    y[[dist.name]] <- as.data.frame(y[[dist.name]])
    y[[dist.name]]$type <- grp
    centroids <- aggregate(cbind(PC1,PC2)~type,data=y[[dist.name]],mean)
    y[[dist.name]] <- merge(y[[dist.name]], centroids, by="type", suffixes=c("",".centroid"))
  }
  mtest <- melt(y, id=c('PC1', 'PC2', 'PC1.centroid', 'PC2.centroid', measure=c('type')))
  test <- y
  obj <- ggplot(mtest, aes(x=PC1, y=PC2, color=type)) + 
    geom_point(size=3) +
    geom_point(data=centroids,aes(x=PC1,y=PC2,color=type)) +
    geom_segment(aes(x=PC1.centroid, y=PC2.centroid, xend=PC1, yend=PC2, color=type)) +
    stat_ellipse() +
    facet_wrap(~ L1, scale="free") + theme_bw()
  return(obj)
}

generate_distance_barplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                       grp.name, strata=NULL, within=T, between=T, bt.no = 100, ann='') {
  obj1 <- NULL
  strata.name <- strata
  df <- data.obj$meta.dat
  grp <- as.factor(df[, grp.name])
  grp.levels <- levels(grp)
  grp.nlevels <- nlevels(grp)
  grp.btws <- outer(grp.levels, grp.levels, paste, sep="#")
  grp.btws <- grp.btws[lower.tri(grp.btws)]
  
  if (!is.null(strata.name)) {
    strata <- df[, strata.name]
  } else {
    strata <- factor(rep(1, nrow(df))) # pseudo strata
  }
  res.df <- NULL
  for (dist.name in dist.names) {
    for (stratum in levels(strata)) {
      ind <- strata %in% stratum
      dist.sub <- dist.obj[[dist.name]][ind, ind]
      df2 <- df[ind, , drop=FALSE]
      if (between) {
        for (grp.btw in grp.btws) {
          ind1 <- which(df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[1])
          ind2 <- which(df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[2])	
          temp <- as.vector(dist.sub[ind1, ind2])
          sem <- sd(sapply(1:bt.no, function (i) {
            mean(dist.sub[as.numeric(sample(paste(ind1), repl=TRUE)), 
                          as.numeric(sample(paste(ind2), repl=TRUE))])
          }))
          res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Between', 
                                             DistanceType=grp.btw, Distance=mean(temp), sd=sem))
        }	
      }
      
      if (within) {
        for (grp.wth in grp.levels) {
          ind1 <- which(df2[, grp.name] == grp.wth)
          temp <- dist.sub[ind1, ind1]			
          temp <- temp[lower.tri(temp)]
          sem <- sd(sapply(1:bt.no, function (i) {
            ind2 <- as.numeric(sample(paste(ind1), repl = TRUE))
            temp <- dist.sub[ind2, ind2]			
            temp <- temp[lower.tri(temp)]
            mean(temp)
          }))
          res.df <- rbind(res.df, data.frame(DistanceMetric=dist.name, Strata=stratum, Compare='Within',
                                             DistanceType=grp.wth, Distance=mean(temp), sd=sem))
        }	
      }
    }
    
  }
  if (between & within) {
    res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
    levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
  } else {
    if (between) {
      res.df$DistanceType <- factor(res.df$DistanceType)
      levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'))
    } else {
      res.df$DistanceType <- factor(res.df$DistanceType)
      levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
    }
  }
  
  if (is.null(strata.name)) {
    limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
    dodge <- position_dodge(width=0.9)
    for (dist.name in dist.names) {
      cat(dist.name, "...\n")
      temp <- res.df[res.df$DistanceMetric == dist.name, ]
      obj1 <- ggplot(temp, aes(x=DistanceType, y=Distance, fill=DistanceType)) + 
        geom_bar(position=dodge, stat="identity", width=0.75) + 
        geom_bar(position=dodge, stat="identity", width=0.75, colour="black", show_guide=FALSE, size=0.25) +
        geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
        labs(y=paste(dist.name, "Distance"), x='') +
        theme(legend.position="none") +
        theme(axis.text.x=element_text(angle=90, hjust=1))
      
    }
  } else {
    limits <- aes(ymax = Distance + sd, ymin=Distance - sd)
    dodge <- position_dodge(width=0.9)
    for (dist.name in dist.names) {
      cat(dist.name, "...\n")
      temp <- res.df[res.df$DistanceMetric == dist.name, ]
      obj1 <- ggplot(temp, aes(x=Strata, y=Distance, fill=DistanceType)) + 
        geom_bar(position=dodge, stat="identity") + 
        geom_bar(position=dodge, stat="identity", colour="black", show_guide=FALSE, size=0.25) +
        geom_errorbar(limits, position=dodge, size=0.25, width=0.5) +
        labs(y=paste(dist.name, "Distance"), x=strata.name) +
        theme(axis.text.x=element_text(angle=90, hjust=1))
    }
  }
  return(obj1)
}

generate_distance_boxplot <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                       grp.name, strata=NULL, within=F, between=T, ann='') {
  obj1 <- NULL
  strata.name <- strata
  df <- data.obj$meta.dat
  grp <- as.factor(df[, grp.name])
  grp.levels <- levels(grp)
  grp.nlevels <- nlevels(grp)
  grp.btws <- outer(grp.levels, grp.levels, paste, sep="#")
  grp.btws <- grp.btws[lower.tri(grp.btws)]
  
  if (!is.null(strata.name)) {
    strata <- df[, strata.name]
  } else {
    strata <- factor(rep(1, nrow(df))) # pseudo strata
  }
  res.df <- NULL
  for (dist.name in dist.names) {
    for (stratum in levels(strata)) {
      ind <- strata %in% stratum
      dist.sub <- dist.obj[[dist.name]][ind, ind]
      df2 <- df[ind, , drop=FALSE]
      if (between) {
        for (grp.btw in grp.btws) {
          ind1 <- df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[1]
          ind2 <- df2[, grp.name] == unlist(strsplit(grp.btw, "#"))[2]	
          temp <- as.vector(dist.sub[ind1, ind2])
          n <- length(temp)
          res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Between', n),
                                             DistanceType=rep(grp.btw, n), Distance=temp))
        }	
      }
      if (within) {
        for (grp.wth in grp.levels) {
          ind1 <- df2[, grp.name] == grp.wth
          temp <- dist.sub[ind1, ind1]			
          temp <- temp[lower.tri(temp)]
          n <- length(temp)
          res.df <- rbind(res.df, data.frame(DistanceMetric=rep(dist.name, n), Strata=rep(stratum, n), Compare=rep('Within', n),
                                             DistanceType=rep(grp.wth, n), Distance=temp))
        }	
      }
    }
    
  }
  if (between & within) {
    res.df$DistanceType <- factor(res.df$DistanceType, levels=c(grp.btws, grp.levels))
    levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'), paste('Within\n', grp.levels))
  } else {
    if (between) {
      res.df$DistanceType <- factor(res.df$DistanceType)
      levels(res.df$DistanceType) <- c(sapply(strsplit(grp.btws, "#"), paste, collapse=' vs\n'))
    } else {
      res.df$DistanceType <- factor(res.df$DistanceType)
      levels(res.df$DistanceType) <- c(paste('Within\n', grp.levels))
    }
  }
  
  if (is.null(strata.name)) {
    dodge <- position_dodge(width=0.95)		
    obj1 <- ggplot(res.df, aes(x=DistanceType, y=Distance, col=DistanceType)) + 
      geom_boxplot(position=dodge, outlier.colour = NA) + 
      #					geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1)) +
      labs(x='') +
      theme(legend.position="none") + facet_wrap(~DistanceMetric)
    #}
  } else {
    dodge <- position_dodge(width=0.95)		
    obj1 <- ggplot(res.df, aes(x=Strata, y=Distance, col=DistanceType)) + 
      geom_boxplot(position=dodge, outlier.colour = NA) + 
      labs(x=strata.name) + facet_wrap(~DistanceMetric)
  }
  return(obj1)
}

# Rev: 2016_09_10, Implement p value based omnibus test
PermanovaG2 <- function (formula, dat = NULL, ...) 
{
  save.seed <- get(".Random.seed", .GlobalEnv)
  lhs <- formula[[2]]
  lhs <- eval(lhs, dat, parent.frame())
  rhs <- as.character(formula)[3]
  p.perms <- list()
  p.obs <- list()
  for (i in 1:(dim(lhs)[3])) {
    assign(".Random.seed", save.seed, .GlobalEnv)
    Y <- as.dist(lhs[, , i])
    formula2 <- as.formula(paste("Y", "~", rhs))
    obj <- adonis(formula2, dat, ...)
    perm.mat <- obj$f.perms
    p.perms[[i]] <- 1 - (apply(perm.mat, 2, rank) - 1) / nrow(perm.mat)
    p.obs[[i]] <- obj$aov.tab[1:ncol(perm.mat), "Pr(>F)"]
    
  }
  
  omni.pv <- NULL
  indiv.pv <- NULL
  for (j in 1:ncol(perm.mat)) {
    p.perms.j <- sapply(p.perms, function (x) x[, j])
    p.obj.j <- sapply(p.obs, function (x) x[j])
    omni.pv <- c(omni.pv, mean(c(rowMins(p.perms.j ) <= min(p.obj.j), 1)))
    indiv.pv <- rbind(indiv.pv, p.obj.j)
  }
  colnames(indiv.pv) <- paste0('D', 1:ncol(indiv.pv), '.p.value')
  rownames(indiv.pv) <- 1:nrow(indiv.pv)
  
  aov.tab <- data.frame(indiv.pv, omni.p.value = omni.pv)
  rownames(aov.tab) <- rownames(obj$aov.tab)[1:ncol(perm.mat)]
  list(aov.tab = aov.tab)
}

perform_permanova_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
                                    PermanovaG.dist=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'),
                                    formula=NULL,  grp.name=NULL, adj.name=NULL, pairwise=F, block.perm=F, strata=NULL, ann='', ...) {
  # PermanovaG not implemented for block permutation
  result <- list()
  
  df <- data.obj$meta.dat
  if (!is.null(strata)) {
    if (is.character(strata)) {
      strata <- df[, strata]
    }
  }
  if (!is.null(formula)) {
    ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
    df <- df[ind, ]
    if (!is.null(strata)) {
      strata <- strata[ind]
    }

    permanova.obj <- list()
    for (dist.name in dist.names) {

      dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
      if (block.perm == F) {
        obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
      } else {
        obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
      }
      
      prmatrix(obj$aov.tab)
      permanova.obj[[dist.name]] <- obj$aov.tab

    }
    result$permanova.obj <- permanova.obj
    permanovaG.obj <- NULL
    if (block.perm == F & !is.null(PermanovaG.dist)) {
      response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
      for (dist.name in PermanovaG.dist) {
        response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
      }
      obj <- PermanovaG2(as.formula(paste("response", formula)), df,  strata=strata, ...)
      prmatrix(obj$aov.tab)
      permanovaG.obj <- obj$aov.tab

      result$permanovaG.obj <- permanovaG.obj
    }
    cat("\n")
    sink()
  } else {
    if (pairwise == F) {
      if (is.null(adj.name)) {
        formula <- paste('~', grp.name)
      } else {
        formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
      }
      
      ind <- apply(df[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
      df <- df[ind, ]
      if (!is.null(strata)) {
        strata <- strata[ind]
      }

      permanova.obj <- list()
      for (dist.name in dist.names) {

        dist.mat <- as.dist(dist.obj[[dist.name]][ind, ind])
        if (block.perm == F) {
          obj <- adonis(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
        } else {
          obj <- adonis2(as.formula(paste("dist.mat", formula)), df,  strata=strata, ...)
        }
        prmatrix(obj$aov.tab)
        permanova.obj[[dist.name]] <- obj$aov.tab

      }
      result$permanova.obj <- permanova.obj
      permanovaG.obj <- NULL
      if (block.perm == F & !is.null(PermanovaG.dist)) {
        response <- array(NA, c(sum(ind), sum(ind), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
        for (dist.name in PermanovaG.dist) {
          response[, , dist.name] <- dist.obj[[dist.name]][ind, ind]
        }
        obj <- PermanovaG2(as.formula(paste("response", formula)), df,  strata=strata, ...)
        prmatrix(obj$aov.tab)
        permanovaG.obj <- obj$aov.tab
        cat("\n")
        result$permanovaG.obj <- permanovaG.obj
      }
      
    } else {
      
      grp <- factor(df[, grp.name])
      grp.levels <- levels(grp)
      grp.nlevels <- nlevels(grp)
      pmat.all <- NULL
      rmat.all <- NULL
      pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
      colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
      for (dist.name in dist.names) {
        cat(dist.name, " distance: \n")
        pmat <- matrix(NA, grp.nlevels, grp.nlevels)
        colnames(pmat) <- rownames(pmat) <- grp.levels
        rmat <- matrix(NA, grp.nlevels, grp.nlevels)
        colnames(rmat) <- rownames(rmat) <- grp.levels
        for (i in 1:(grp.nlevels-1)) {
          grp.level1 <- grp.levels[i]
          for (j in (i+1):grp.nlevels) {
            
            grp.level2 <- grp.levels[j]
            cat(grp.level1, ' vs ', grp.level2, '\n')
            ind <- grp %in% c(grp.level1, grp.level2)
            df2 <- subset(df, ind)
            df2[, grp.name] <- factor(df2[, grp.name])
            dist.mat <- dist.obj[[dist.name]][ind, ind]
            strata2 <- strata[ind]
            
            if (is.null(adj.name)) {
              formula <- paste('~', grp.name)
            } else {
              formula <- 	paste('~', paste(adj.name, collapse='+'), '+', grp.name)		
            }
            
            ind2 <- apply(df2[, gsub("^\\s+|\\s+$", "", strsplit(strsplit(formula, "~")[[1]][2], "\\+")[[1]]), drop=F], 1, function(x) sum(is.na(x))) == 0
            df2 <- df2[ind2, ]
            dist.mat2 <- as.dist(dist.mat[ind2, ind2])
            strata2 <- strata2[ind2]
            if (block.perm == F) {
              obj <- adonis(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2, ...)
            } else {
              obj <- adonis2(as.formula(paste("dist.mat2", formula)), df2,  strata=strata2,  ...)
            }
            prmatrix(obj$aov.tab)
            cat("\n")
            
            if (block.perm == F) {
              pmat[i, j] <- pmat[j, i] <- obj$aov.tab[length(adj.name)+1, 6]
              rmat[i, j] <- rmat[j, i] <- obj$aov.tab[length(adj.name)+1, 5]
            } else {
              pmat[i, j] <- pmat[j, i] <- obj$aov.tab[1, 6]
              rmat[i, j] <- rmat[j, i] <- obj$aov.tab[1, 5]
            }
            
            # PERMANOVA G after last distance
            if (block.perm == F & !is.null(PermanovaG.dist)) {
              if (dist.name == dist.names[length(dist.names)]) {
                cat('\nPERMANOVA G test combining ', paste(PermanovaG.dist, collapse=','), '\n')
                response <- array(NA, c(sum(ind2), sum(ind2), length(PermanovaG.dist)), dimnames=list(NULL, NULL, PermanovaG.dist))
                for (dist.name in PermanovaG.dist) {
                  response[, , dist.name] <- dist.mat[ind2, ind2]
                }
                obj <- PermanovaG2(as.formula(paste("response", formula)), df2,  strata=strata2, ...)
                prmatrix(obj$aov.tab)

                pmat.G[i, j] <- pmat.G[j, i] <- obj$aov.tab[length(adj.name)+1, 'omni.p.value']
              }
            }
          }
        }

        pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
        rmat.all <- rbind(rmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(rmat), rep("", grp.nlevels))
      }
      
      result$pmat.all <- pmat.all
      result$rmat.all <- rmat.all
      result$pmat.G <- pmat.G
      
    }
  }
  
  return(result)
}

# Rev: 2016_12_02, MiKRAT for binary result, add out_type='D'
perform_mirkat_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), 
                                 grp.name=NULL, adj.name=NULL, pairwise=F,  ann='', ...) {
  
  # MiRKAT not implemented for correlated data
  df <- data.obj$meta.dat
  result <- list()
  if (pairwise == F) {
    
    ind <- apply(df[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
    df <- df[ind, ]
    
    grp <- df[, grp.name]
    
    if (is.character(grp)) {
      grp <- factor(grp)
    }
    # Rev: 2016_12_02
    if (is.factor(grp)) {
      if (nlevels(grp) > 2) {
        stop('Currently MiRKAT only supports binary outcome!')
      } else {
        grp <- as.numeric(grp) - 1
        out_type <- 'D'
      }
    } else {
      out_type <- 'C'
    }
    if (!is.null(adj.name)) {
      # No intercept
      adj <- model.matrix(~ ., data.frame(df[, adj.name]))
      # Remove collinear terms
      qadj <- qr(adj, tol = 1e-07)
      adj <- adj[, qadj$pivot, drop = FALSE]
      adj <- adj[, 1:qadj$rank, drop = FALSE]
      # Remove intercept
      adj <- adj[, colSums(adj==1) != nrow(adj)]
    } else {
      adj <- NULL
    }
    
   
    Ks <- list()
    for (dist.name in dist.names) {
      Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][ind, ind])
    }
    # Rev: 2016_12_02
    obj <- MiRKAT(grp, X=adj, Ks, out_type=out_type)
    
    result$indiv <- prmatrix(t(obj$indivP))
    result$omni <- obj$omnibus_p

  } else {
    
    grp <- factor(df[, grp.name])
    grp.levels <- levels(grp)
    grp.nlevels <- nlevels(grp)
    pmat.G <- matrix(NA, grp.nlevels, grp.nlevels)
    colnames(pmat.G) <- rownames(pmat.G) <- grp.levels
    parr <- array(NA, c(grp.nlevels, grp.nlevels, length(dist.names)), dimnames=list(grp.levels, grp.levels, dist.names))
    
    for (i in 1:(grp.nlevels-1)) {
      grp.level1 <- grp.levels[i]
      for (j in (i+1):grp.nlevels) {
        
        grp.level2 <- grp.levels[j]
        cat(grp.level1, ' vs ', grp.level2, '\n')
        ind <- grp %in% c(grp.level1, grp.level2)
        df2 <- subset(df, ind)
        df2[, grp.name] <- factor(df2[, grp.name])
        ind2 <- apply(df2[, c(grp.name, adj.name), drop=F], 1, function(x) sum(is.na(x))) == 0
        df2 <- df[ind2, ]
        
        grp2 <- as.numeric(df2[, grp.name]) - 1
        if (!is.null(adj.name)) {
          # No intercept
          adj <- model.matrix(~ ., data.frame(df2[, adj.name]))
          # Remove collinear terms
          qadj <- qr(adj, tol = 1e-07)
          adj <- adj[, qadj$pivot, drop = FALSE]
          adj <- adj[, 1:qadj$rank, drop = FALSE]
          # Remove intercept
          adj <- adj[, colSums(adj==1) != nrow(adj)]
        } else {
          adj <- NULL
        }
        
        Ks <- list()
        for (dist.name in dist.names) {
          Ks[[dist.name]] <- D2K(dist.obj[[dist.name]][rownames(df2), rownames(df2)])
        }
        # Rev: 2016_12_02
        obj <- MiRKAT(grp2, X=adj, Ks, out_type='D')
        pmat.G[i, j] <- pmat.G[j, i] <- obj$omnibus_p
        parr[i, j, ] <- parr[j, i, ] <- obj$indivP

        result$indiv <- prmatrix(t(obj$indivP))

        result$omni <- obj$omnibus_p

        
      }
    }
    pmat.all <- NULL
    for (dist.name in dist.names) {
      pmat <- parr[, , dist.name]
      pmat.all <- rbind(pmat.all, c(dist.name, rep('', grp.nlevels-1)), formatC(pmat), rep("", grp.nlevels))
    }
    
    result$pmat.all <- pmat.all
    result$pmat.G <- pmat.G
  }
  return(result)
}

perform_betadisper_test <- function (data.obj, dist.obj, dist.names=c('UniFrac', 'GUniFrac', 'WUniFrac', 'BC'), grp.name) {
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  result <- list()
  for (dist.name in dist.names) {

    dist.mat <- as.dist(dist.obj[[dist.name]])
    obj <- betadisper(dist.mat, grp)
    result[[dist.name]] <- prmatrix(anova(obj))
  }
  return(result)
}

twopart.test <- function(x1, x2, zero.p=0.2) {
  
  n1 <- length(x1)
  n2 <- length(x2)
  p1 <- mean(x1 != 0)
  p2 <- mean(x2 != 0)
  m1 <- sum(x1 != 0)
  m2 <- sum(x2 != 0)
  p12 <- (m1 + m2) / (n1 + n2)
  q12 <- 1 - p12
  
  if (q12 >= zero.p) {
    Z <- (abs(p1 - p2) - (1/(2*n1) + 1/(2*n2))) / sqrt(p12 * q12 * (1/n1 + 1/n2))
    x1 <- x1[x1!=0]
    x2 <- x2[x2!=0]
    R1 <- sum(rank(c(x1, x2))[1:length(x1)])
    ti <- as.vector(table(c(x1, x2)))
    W <- (abs(R1 - m1*(m1+m2+1)/2) - 1/2) / sqrt((m1*m2/12)*(m1+m2+1-sum(ti*(ti^2-1))/(m1+m2)/(m1+m2-1)))
    X2 <- Z^2 + W^2
    res <- list()
    res$stat <- X2
    res$p.value <- 1 - pchisq(X2, 2)
    res$Z <- Z
    res$W <- W
    res$test <- 'TwoPart'
  } else {
    res <- wilcox.test(x1, x2)
    res$test <- 'Wilcox'
  }
  res
}

getPermuteMatrix <- function (perm, N, strata = NULL) 
{
  if (length(perm) == 1) {
    perm <- how(nperm = perm)
  }
  if (!missing(strata) && !is.null(strata)) {
    if (inherits(perm, "how") && is.null(getBlocks(perm))) 
      setBlocks(perm) <- strata
  }
  if (inherits(perm, "how")) 
    perm <- shuffleSet(N, control = perm)
  if (is.null(attr(perm, "control"))) 
    attr(perm, "control") <- structure(list(within = list(type = "supplied matrix"), 
                                            nperm = nrow(perm)), class = "how")
  perm
}

perm_fdr_adj <- function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  perm.no <- ncol(Fp)
  Fp <- as.vector(Fp)
  Fp <- Fp[!is.na(Fp)]
  Fp <- sort(c(Fp, F0), decreasing = F)
  n <- length(Fp)
  m <- length(F0)
  FPN <- (n + 1) - match(F0, Fp) - 1:m
  p.adj.fdr <- FPN / perm.no / (1:m)
  #		p.adj.fdr <- sapply(F0, function(x) sum(Fp >= 
  #									x, na.rm=TRUE) / perm.no)/(1:length(F0))
  p.adj.fdr <- pmin(1, rev(cummin(rev(p.adj.fdr))))[order(ord)]
}

perm_fwer_adj <- function (F0, Fp) {
  ord <- order(F0, decreasing = T)
  F0 <- F0[ord]
  col.max <- colMaxs(Fp, na.rm=TRUE)
  p.adj.fwer <- sapply(F0, function(x) mean(col.max >= x))[order(ord)]
}

# Need to be further comprehensively tested
# a. Speed up, finished!
# b. Permutation method （response, covariate or residual permutation), residual will be default!
# c. Revise - add permutation stratified by subject, finished!
# Rev: 2017_02_02 permutation-based FDR control
# Rev: 2017_02_13 Add LMM-based permutation for block.perm=TRUE (type I error is controled, but power study hasn't been comprehensively studied)
# Rev: 2017_02_24 allow 'adj.name' to contain multiple covariates
# Still need to address NA's, currently simply remove NA's
permute_differential_analysis <- function (meta.dat, prop, grp.name, adj.name=NULL, strata=NULL, 
                                           block.perm=FALSE, sqrt.trans=TRUE, resid.perm=TRUE, perm.no=999) {
  # Square root transformation
  # User should take care of the normalization, transformation and addressing outliers
  if (sqrt.trans) {
    Y <- sqrt(prop)
  } else {
    Y <- prop
  }
  row.names <- rownames(Y)
  
  if (!is.null(strata)) {
    strata <- factor(strata)
  }
  
  # Prepare model matrix
  n <- ncol(prop)
  I <- diag(n)
  if (is.null(adj.name)) {
    M0 <- model.matrix(~ 1, meta.dat)
  } else {
    df0 <- meta.dat[, c(adj.name), drop=F]
    M0 <- model.matrix( ~., df0)
  }
  
  #	P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
  
  df1 <- meta.dat[, c(adj.name, grp.name), drop=F]
  M1 <- model.matrix( ~., df1)
  
  # QR decompostion
  qrX0 <- qr(M0, tol = 1e-07)
  Q0 <- qr.Q(qrX0)
  Q0 <- Q0[, 1:qrX0$rank, drop=FALSE]
  
  qrX1 <- qr(M1, tol = 1e-07)
  Q1 <- qr.Q(qrX1)
  Q1 <- Q1[, 1:qrX1$rank, drop=FALSE]
  
  # Got residual
  if (resid.perm) {
    if (!block.perm) {
      # Permute the residual
      if (is.null(adj.name)) {
        Y <- t(resid(lm(as.formula(paste('t(Y) ~ 1')), meta.dat)))
      } else {
        Y <- t(resid(lm(as.formula(paste('t(Y) ~ ', paste(adj.name, collapse='+'))), meta.dat)))
      }
      
    } else {
      
      if (is.null(strata)) {
        stop('Block permutation requires strata!\n')
      } else {
        Y.r <- matrix(NA, nrow(Y), nlevels(strata))
        Y.e <- Y
        cat('Fitting linear mixed effects model ...\n')
        for (j in 1:nrow(Y)) {
          # Linear mixed effects model
          yy <- Y[j, ]
          meta.dat$yy <- yy
          meta.dat$strata <- strata
          if (is.null(adj.name)) {
            obj <- lme(as.formula(paste('yy ~ 1')), random =~ 1 |  strata, data=meta.dat, method='ML')
            # The order is the same as the levels
            Y.r[j, ] <- random.effects(obj)[, 1]
            Y.e[j, ] <- resid(obj)
          } else {
            obj <- lme(as.formula(paste('yy ~ ', paste(adj.name, collapse='+'))), random =~ 1 |  strata, data=meta.dat, method='ML')
            Y.r[j, ] <- random.effects(obj)[, 1]
            Y.e[j, ] <- resid(obj)
          }
        }
        Y <- Y.r[, as.numeric(strata)] + Y.e
        #			Y <- Y - rowMeans(Y)
        #		    Y <- Y.e
      }
    }
  }
  
  
  
  TSS <- rowSums(Y^2)
  MSS1 <- rowSums((Y %*% Q1)^2)
  MSS0 <- rowSums((Y %*% Q0)^2)  # Not necessary, it's zero
  F0 <- (MSS1 - MSS0) /  (TSS - MSS1) 
  
  #	P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
  #	F0 <- diag(Y %*% (P1 - P0) %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y))
  #	df3 <- df1
  
  if (block.perm == FALSE) {
    perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
    perm.no <- nrow(perm.ind)
  }
  
  cat('Permutation test ....\n')
  Fp <- sapply(1:perm.no, function(i) {
    if (i %% 100 == 0) cat('.')
    if (block.perm == FALSE) {
      Yp <- Y[, perm.ind[i, ]]
    } else {
      # Double permutation
      strata.p <- factor(strata, levels=sample(levels(strata)))
      Yp <- Y.r[, as.numeric(strata.p)] + Y.e[, sample(ncol(Y))]
      #				    Yp <- Y.e[, sample(ncol(Y))]
      #					Yp <- Yp - rowMeans(Yp)
    }
    
    #				df3[, grp.name] <- sample(df1[, grp.name])
    #				M1 <- model.matrix( ~., df3)
    #				P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
    MSS1p <- rowSums((Yp %*% Q1)^2)
    MSS0p <- rowSums((Yp %*% Q0)^2)    
    if (block.perm == FALSE) {
      TSSp <- TSS
    } else {
      TSSp <- rowSums(Yp^2)
    }
    (MSS1p - MSS0p) /  (TSSp - MSS1p) 
  })
  
  
  if (mean(is.na(F0)) >= 0.1) {
    warning('More than 10% observed F stats have NA! Please check! \n')
  }
  
  if (mean(is.na(Fp)) >= 0.1) {
    warning('More than 10% permuted F stats have NA! Please check! \n')
  }
  
  na.ind <- is.na(F0)
  F0 <- F0[!na.ind]
  Fp <- Fp[!na.ind, ]
  
  p.raw <- cbind(Fp >= F0, 1)
  p.raw <- rowMeans(p.raw)
  #	p.raw[is.na(p.raw)] <- 1   
  
  p.adj.fdr <- perm_fdr_adj(F0, Fp)
  p.adj.fwer <- perm_fwer_adj(F0, Fp)
  
  # Pad back the NA values
  pad <- function (vec, ind) {
    vec0 <- numeric(length(ind))
    vec0[!ind] <- vec
    vec0[ind] <- NA
    vec0
  }
  
  F0 <- pad(F0, na.ind)
  p.raw <- pad(p.raw, na.ind)
  p.adj.fdr <- pad(p.adj.fdr, na.ind)
  p.adj.fwer <- pad(p.adj.fwer, na.ind)
  
  names(F0) <- names(p.raw) <- names(p.adj.fdr) <- names(p.adj.fwer) <- row.names
  return(list(F.stat=F0, p.raw=p.raw, p.adj.fdr=p.adj.fdr, p.adj.fwer=p.adj.fwer))
  #	return(p.raw)
}
#


# New: 2017_08_17 Add a new variant of PERMANOVA with matrix decomposition
# Partially validated
permanova2 <- function (meta.dat, D, grp.name, adj.name=NULL, strata=NULL, 
                        block.perm=FALSE, resid.perm=TRUE, perm.no=999, eig = c('All', 'Positive')) {
  
  eig <- match.arg(eig)
  # Square root transformation
  # User should take care of the normalization, transformation and addressing outliers
  
  D <- as.matrix(D)^2
  n <- nrow(D)
  G <- mean(D) + D - rowMeans(D) - matrix(rep(1, n), ncol=1) %*% colMeans(D)
  
  eig.obj <- eigen(G, symm=TRUE)
  
  if (eig == 'All') {
    Y <- eig.obj$vectors
    lambda <- eig.obj$values
  }
  
  if (eig == 'Positive') {
    lambda <- eig.obj$values
    ind <- lambda > 1e-6
    lambda <- lambda[ind]
    Y <- Y[, ind, drop = FALSE]
  }
  
  Y <- t(Y)
  
  if (!is.null(strata)) {
    strata <- factor(strata)
  }
  
  # Prepare model matrix
  n <- ncol(Y)
  I <- diag(n)
  if (is.null(adj.name)) {
    M0 <- model.matrix(~ 1, meta.dat)
  } else {
    df0 <- meta.dat[, c(adj.name), drop=F]
    M0 <- model.matrix( ~., df0)
  }
  
  #	P0 <- M0 %*% solve(t(M0) %*% M0) %*% t(M0)
  
  df1 <- meta.dat[, c(adj.name, grp.name), drop=F]
  M1 <- model.matrix( ~., df1)
  
  # QR decompostion
  qrX0 <- qr(M0, tol = 1e-07)
  Q0 <- qr.Q(qrX0)
  Q0 <- Q0[, 1:qrX0$rank, drop=FALSE]
  
  qrX1 <- qr(M1, tol = 1e-07)
  Q1 <- qr.Q(qrX1)
  Q1 <- Q1[, 1:qrX1$rank, drop=FALSE]
  
  # Got residual
  if (resid.perm) {
    if (!block.perm) {
      # Permute the residual
      if (is.null(adj.name)) {
        Y <- t(resid(lm(as.formula(paste('t(Y) ~ 1')), meta.dat)))
      } else {
        Y <- t(resid(lm(as.formula(paste('t(Y) ~ ', paste(adj.name, collapse='+'))), meta.dat)))
      }
      
    } else {
      
      if (is.null(strata)) {
        stop('Block permutation requires strata!\n')
      } else {
        Y.r <- matrix(NA, nrow(Y), nlevels(strata))
        Y.e <- Y
        #				cat('Fitting linear mixed effects model ...\n')
        for (j in 1:nrow(Y)) {
          # Linear mixed effects model
          yy <- Y[j, ]
          meta.dat$yy <- yy
          meta.dat$strata <- strata
          if (is.null(adj.name)) {
            obj <- lme(as.formula(paste('yy ~ 1')), random =~ 1 |  strata, data=meta.dat, method='ML')
            # The order is the same as the levels
            Y.r[j, ] <- random.effects(obj)[, 1]
            Y.e[j, ] <- resid(obj)
          } else {
            obj <- lme(as.formula(paste('yy ~ ', paste(adj.name, collapse='+'))), random =~ 1 |  strata, data=meta.dat, method='ML')
            Y.r[j, ] <- random.effects(obj)[, 1]
            Y.e[j, ] <- resid(obj)
          }
        }
        Y <- Y.r[, as.numeric(strata)] + Y.e
        #			Y <- Y - rowMeans(Y)
        #		    Y <- Y.e
      }
    }
  }
  
  
  
  TSS <- rowSums(Y^2)
  MSS1 <- rowSums((Y %*% Q1)^2)
  MSS0 <- rowSums((Y %*% Q0)^2)  # Not necessary, it's zero
  F0 <- sum(lambda * (MSS1 - MSS0)) /  sum(lambda * (TSS - MSS1))
  
  #	P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
  #	F0 <- diag(Y %*% (P1 - P0) %*% t(Y)) / diag(Y %*% (I - P1) %*% t(Y))
  #	df3 <- df1
  
  if (block.perm == FALSE) {
    perm.ind <- vegan:::getPermuteMatrix(perm.no, n, strata = strata)
    perm.no <- nrow(perm.ind)
  }
  
  #	cat('Permutation test ....\n')
  Fp <- sapply(1:perm.no, function(i) {
    #				if (i %% 100 == 0) cat('.')
    if (block.perm == FALSE) {
      Yp <- Y[, perm.ind[i, ]]
    } else {
      # Double permutation
      strata.p <- factor(strata, levels=sample(levels(strata)))
      Yp <- Y.r[, as.numeric(strata.p)] + Y.e[, sample(ncol(Y))]
      #				    Yp <- Y.e[, sample(ncol(Y))]
      #					Yp <- Yp - rowMeans(Yp)
    }
    
    #				df3[, grp.name] <- sample(df1[, grp.name])
    #				M1 <- model.matrix( ~., df3)
    #				P1 <- M1 %*% solve(t(M1) %*% M1) %*% t(M1)
    MSS1p <- rowSums((Yp %*% Q1)^2)
    MSS0p <- rowSums((Yp %*% Q0)^2)    
    if (block.perm == FALSE) {
      TSSp <- TSS
    } else {
      TSSp <- rowSums(Yp^2)
    }
    sum(lambda * (MSS1p - MSS0p)) / sum(lambda * (TSSp - MSS1p)) 
  })
  
  p.value <- mean(c(Fp >= F0, 1))
  
  return(list(f0 = F0, f.perms = Fp, p.value = p.value))
  #	return(p.raw)
}


# New: 2017_02_07
GMPR <- function (comm, intersect.no=4, ct.min=2) {
  # Computes the GMPR size factor
  #
  # Args:
  #   comm: a matrix of counts, row - features (OTUs, genes, etc) , column - sample
  #   intersect.no: the minimum number of shared features between sample pair, where the ratio is calculated
  #   ct.min: the minimum number of counts required to calculate ratios ct.min = 5 has better results
  
  #
  # Returns:
  #   a list that contains:
  #      gmpr： the GMPR size factors for all samples; Samples with distinct sets of features will be output as NA.
  #      nss:   number of samples with significant sharing (> intersect.no) including itself
  
  # mask counts < ct.min
  comm[comm < ct.min] <- 0
  
  if (is.null(colnames(comm))) {
    colnames(comm) <- paste0('S', 1:ncol(comm))
  }
  
  cat('Begin GMPR size factor calculation ...\n')
  
  comm.no <- numeric(ncol(comm))
  gmpr <- sapply(1:ncol(comm),  function(i) {		
    if (i %% 50 == 0) {
      cat(i, '\n')
    }
    x <- comm[, i]
    # Compute the pairwise ratio
    pr <- x / comm
    # Handling of the NA, NaN, Inf
    pr[is.nan(pr) | !is.finite(pr) | pr == 0] <- NA
    # Counting the number of non-NA, NaN, Inf
    incl.no <- colSums(!is.na(pr))		
    # Calculate the median of PR
    pr.median <- colMedians(pr, na.rm=TRUE)
    # Record the number of samples used for calculating the GMPR
    comm.no[i] <<- sum(incl.no >= intersect.no)
    # Geometric mean of PR median
    if (comm.no[i] > 1) {
      return(exp(mean(log(pr.median[incl.no >= intersect.no]))))
    } else {
      return(NA)
    }
  }
  )
  
  if (sum(is.na(gmpr))) {
    warning(paste0('The following samples\n ', paste(colnames(comm)[is.na(gmpr)], collapse='\n'), 
                   '\ndo not share at least ', intersect.no, ' common taxa with the rest samples! ',
                   'For these samples, their size factors are set to be NA! \n', 
                   'You may consider removing these samples since they are potentially outliers or negative controls!\n',
                   'You may also consider decreasing the minimum number of intersecting taxa and rerun the procedure!\n'))
  }
  
  cat('Completed!\n')
  cat('Please watch for the samples with limited sharing with other samples based on NSS! They may be outliers! \n')
  attr(gmpr, 'NSS') <- comm.no
  # Rev: 2017_09_07
  gmpr <- gmpr * median(colSums(comm))
  names(gmpr) <- colnames(comm)
  return(gmpr)
}

# This function for nonparametric/permutaiton method
# Rev: 2017_02_16  Add normalization method; Add transformation; Remove rarefaction (only output warnings);
# Rev: 2017_10_30  Support filtering based on coefficient of variation
# Rev: 2018_01_30 Add ct.min and handle NA size factor
perform_differential_analysis <- function (data.obj, grp.name, adj.name=NULL, subject=NULL, 
                                           taxa.levels=c('Phylum', 'Order', 'Class', 'Family', 'Genus', 'Species'),
                                           method='perm', block.perm=FALSE, perm.no=999,
                                           norm='GMPR', norm.level='Species', intersect.no=4, ct.min = 2,
                                           transform='sqrt',
                                           prev=0.1, minp=0.002, medianp=NULL, cv=NULL,
                                           mt.method='fdr', cutoff=0.15, 
                                           ann='', seed=123, ...) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  res.final <- NULL
  ind <- !is.na(grp)
  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  df <- df[ind, ]
  
  if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
    data.obj$abund.list[['Species']] <- data.obj$otu.tab
    rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
                                                         data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
  }
  
  # Test for sequence-depth confounding
  dep <- colSums(data.obj$otu.tab)
  diff.seq.p <- summary(aov(dep ~ grp))[[1]][1, 'Pr(>F)']
  if (!is.na(diff.seq.p) & diff.seq.p <= 0.05) {
    warning(paste0(
      '\nSignificant sequencing depth confounding with the variable of interest!\n',
      'For nonparametric test/permutaiton test, there may be potentially many false postives (for those less prevalent taxa)!\n',
      'Consider performing rarefaction first! (However, rarefaction will not completly solve the problem.)\n',
      'May also try count-based models, which might have better false postive controls!\n'))
  }
  
  # Calculate size.factor
  if (norm == 'Precalculated') {
    size.factor <- data.obj$size.factor
  }
  if (norm == 'GMPR') {
    if (norm.level %in% c('OTU', 'Species')) {
      tab <- data.obj$otu.tab
    } else {
      tab <- data.obj$abund.list[[norm.level]]
    }
    size.factor <- GMPR(tab, intersect.no, ct.min)
    
    # Rev: 2018_01_30
    ind <- !is.na(size.factor)
    data.obj <- subset_data(data.obj, ind)
    grp <- grp[ind]
    df <- df[ind, ]
    size.factor <- size.factor[ind]
  }	
  if (norm == 'TSS') {
    size.factor <- colSums(data.obj$otu.tab)
  }
  
  # Method-dependent processing
  if (is.null(method)) {
    if (nlevels(grp) == 2) {
      method <- 'wilcox'
    } else {
      method <- 'kruskal'
    }
  }
  if (method == 'wilcox.pair') {
    if (nlevels(grp) != 2) stop("Wilcox test requires two groups!\n")
    if (is.null(subject)) stop("Paired wilcox needs subject information!\n")
    subject <- factor(df[, subject])
    ind1 <- ind2 <- NULL
    for(sub in levels(subject)) {
      temp1 <- which(as.numeric(grp) == 1 & subject == sub)
      temp2 <- which(as.numeric(grp) == 2 & subject == sub)
      if (length(temp1) != 0 & length(temp2) != 0) {
        ind1 <- c(ind1, temp1[1])
        ind2 <- c(ind2, temp2[1])
      }
      
    }
  }
  
  if (method == 'perm.pair') {
    if (is.null(subject)) stop("Paired permutation test needs subject information!\n")
    subject <- factor(df[, subject])
  }
  
  if (method == 'perm') {
    if (!is.null(subject)) {
      subject <- factor(df[, subject])
    } 
  }
  
  if (is.factor(grp)) {
    
    pv.list <- qv.list <-  fc.list <- pc.list <- m.list <- nzm.list  <- prv.list <- list()
    res.final <- NULL
    for (LOI in taxa.levels) {
      cat(LOI, "\n")
      ct <- data.obj$abund.list[[LOI]]
      
      # Filtering
      prop0 <- t(t(ct) / colSums(ct))
      if (!is.null(prev)) {
        prop0 <- prop0[rowSums(prop0!=0) > prev * ncol(prop0), , drop=FALSE]	
        ct <- ct[rownames(prop0), , drop=FALSE]
      }
      if (!is.null(minp)) {
        prop0 <- prop0[rowMaxs(prop0) > minp, , drop=FALSE]	
        ct <- ct[rownames(prop0), , drop=FALSE]
      }
      
      if (!is.null(medianp)) {
        nz.mean <- apply(prop0, 1, function(x) median(x[x!=0]))
        prop0 <- prop0[nz.mean > medianp, , drop=FALSE]	
        ct <- ct[rownames(prop0), , drop=FALSE]
      }
      
      if (!is.null(cv)) {
        prop0 <- prop0[rowSds(prop0) / rowMeans(prop0) > cv, , drop=FALSE]	
        ct <- ct[rownames(prop0), , drop=FALSE]
      }
      
      # Normalization
      prop <- t(t(ct) / size.factor)
      
      # Transformation - Others are possible/may be explored in the future
      if (transform == 'sqrt') {
        prop <- sqrt(prop)
      }
      
      if (method == 'perm') {
        set.seed(seed)
        pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, block.perm=block.perm, sqrt.trans=FALSE, perm.no=perm.no)$p.raw
        names(pv.de2) <- rownames(prop)
      }
      
      # For legacy use
      if (method == 'perm.pair') {
        set.seed(seed)
        pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, block.perm=FALSE, sqrt.trans=FALSE, perm.no=perm.no)$p.raw
        names(pv.de2) <- rownames(prop)
      }
      
      pv.vec <- m.vec <- nzm.vec <- prv.vec <-  fc.vec <- pc.vec <-  NULL
      
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        taxon.abund0 <- prop0[taxon, ]
        
        pv <- fc <- pc <- m <- nzm <- prv <- NULL
        if (method == 'wilcox') {
          if (nlevels(grp) != 2) stop("Wilcox test requires two groups!\n")
          pv <- wilcox.test(taxon.abund ~ grp)$p.value
          
        }
        if (method == 'wilcox.pair') {
          pv <- wilcox.test(taxon.abund[ind1], taxon.abund[ind2], paired=T)$p.value
        }
        if (method == 'twopart') {
          if (nlevels(grp) != 2) stop("Two part test requires two groups!\n")
          grp1 <- taxon.abund[as.numeric(grp)==1]
          grp2 <- taxon.abund[as.numeric(grp)==2]
          pv <- twopart.test(grp1, grp2)$p.value
        }	
        if (method == 'kruskal') {
          if (nlevels(grp) <= 2) warning("Kruskal-wallis test requires three or more groups!\n")
          pv <- kruskal.test(taxon.abund ~ grp)$p.value
        }
        
        if (method == 'perm') {
          pv <- pv.de2[taxon]
        }
        
        if (method == 'perm.pair') {
          pv <- pv.de2[taxon]
        }
        m <- tapply(taxon.abund0, grp, function(x) mean(x))
        nzm <- tapply(taxon.abund0, grp, function(x) mean(x[x != 0]))
        prv <- tapply(taxon.abund0, grp, function(x) sum(x != 0))				
        
        # Rev: 2017_02_21 fc change baseline grp
        if (nlevels(grp) == 2) {
          grp.no <- table(grp)
          fc <- log2(m[2] / m[1])
          pc <- prv[2] / grp.no[2] / prv[1] * grp.no[1]
        } else {
          pc <- fc <- NA
        }
        
        pv.vec <- rbind(pv.vec, pv)
        fc.vec <- rbind(fc.vec, fc)
        m.vec <- rbind(m.vec, m)
        nzm.vec <- rbind(nzm.vec, nzm)
        pc.vec <- rbind(pc.vec, pc)
        prv.vec <- rbind(prv.vec, prv / table(grp))	
      }
      
      
      temp <- p.adjust(pv.vec[, 1], 'fdr')
      
      qv.vec <- matrix(temp, ncol=1)
      
      rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(pc.vec) <- rownames(m.vec) <- rownames(nzm.vec) <- rownames(prv.vec) <- rownames(prop)
      colnames(pv.vec) <- 'Pvalue'
      colnames(qv.vec) <- 'Qvalue'
      colnames(fc.vec) <- 'logFoldChange'
      colnames(pc.vec) <- 'PrevalChange'
      colnames(m.vec) <- paste(levels(grp), 'Mean')
      colnames(nzm.vec) <- paste(levels(grp), 'nzMean')
      colnames(prv.vec) <- paste(levels(grp), 'preval')
      
      pv.list[[LOI]] <- pv.vec
      qv.list[[LOI]] <- qv.vec
      fc.list[[LOI]] <- fc.vec
      pc.list[[LOI]] <- pc.vec
      m.list[[LOI]] <- m.vec
      nzm.list[[LOI]] <- nzm.vec
      prv.list[[LOI]] <- prv.vec
      
      res <- cbind(pv.vec, qv.vec, m.vec, nzm.vec, fc.vec, prv.vec, pc.vec)
      rownames(res) <- rownames(prop)
      
      if (mt.method == 'fdr') {
        res.final <- rbind(res.final, res[res[, 'Qvalue'] <= cutoff, , drop=F])
      }
      if (mt.method == 'raw') {
        res.final <- rbind(res.final, res[res[, 'Pvalue'] <= cutoff, , drop=F])
      }
      
    }
    
    if (!is.null(res.final)) {
      colnames(res.final) <- colnames(res)
    }
    return(list(pv.list=pv.list, fc.list=fc.list, pc.list=pc.list, qv.list=qv.list, m.list=m.list, res.final=res.final))
  } else {
    if (is.null(method)) {
      method <- 'Spearman'
    }
    # Continuous case - currently only has DESeq2 
    pv.list <- qv.list <-  fc.list <-  m.list <- list()
    res.final <- NULL
    for (LOI in taxa.levels) {
      cat(LOI, "\n")
      ct <- data.obj$abund.list[[LOI]]
      
      # Filtering
      prop0 <- t(t(ct) / colSums(ct))
      if (!is.null(prev)) {
        prop0 <- prop0[rowSums(prop0!=0) > prev * ncol(prop0), , drop=FALSE]	
        ct <- ct[rownames(prop0), , drop=FALSE]
      }
      if (!is.null(minp)) {
        prop0 <- prop0[rowMaxs(prop0) > minp, , drop=FALSE]	
        ct <- ct[rownames(prop0), , drop=FALSE]
      }
      
      if (!is.null(medianp)) {
        nz.mean <- apply(prop0, 1, function(x) median(x[x!=0]))
        prop0 <- prop0[nz.mean > medianp, , drop=FALSE]	
        ct <- ct[rownames(prop0), , drop=FALSE]
      }
      
      if (!is.null(cv)) {
        prop0 <- prop0[rowSds(prop0) / rowMeans(prop0) > cv, , drop=FALSE]	
        ct <- ct[rownames(prop0), , drop=FALSE]
      }
      
      # Normalization
      prop <- t(t(ct) / size.factor)
      
      # Transformation - Others are possible/may be explored in the future
      if (transform == 'sqrt') {
        prop <- sqrt(prop)
      }
      
      if (method == 'perm') {
        set.seed(seed)
        pv.de2 <- permute_differential_analysis(df, prop, grp.name, adj.name, strata=subject, block.perm=block.perm, sqrt.trans=FALSE, perm.no=perm.no)$p.raw
        names(pv.de2) <- rownames(prop)
        # Place holder
        fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
      }
      
      if (method == 'Spearman') {
        if (!is.null(adj.name)) {
          stop("Spearman test can't adjust covariates!")
        }
        pv.de2 <- apply(prop, 1, function(x) {
          cor.test(x, grp, method='spearman')$p.value
        })
        names(pv.de2) <- rownames(prop)
        # Place holder
        fc.de2 <- sapply(1:nrow(prop), function(i) cor(prop[i, ], df[, grp.name], method='spearman'))
      }
      pv.vec <- matrix(pv.de2, ncol=1)	
      qv.vec <- matrix(p.adjust(pv.vec[, 1], 'fdr'), ncol=1)
      fc.vec <- matrix(fc.de2, ncol=1)
      m.vec <- matrix(rowMeans(prop0), ncol=1)
      
      rownames(pv.vec) <- rownames(qv.vec) <- rownames(fc.vec) <- rownames(m.vec)  <- rownames(prop)
      colnames(pv.vec) <- 'Pvalue'
      colnames(qv.vec) <- 'Qvalue'
      colnames(fc.vec) <- 'SpearmanCorr'
      colnames(m.vec) <- 'Mean'
      
      
      pv.list[[LOI]] <- pv.vec
      qv.list[[LOI]] <- qv.vec
      fc.list[[LOI]] <- fc.vec
      m.list[[LOI]] <- m.vec
      
      
      res <- cbind(m.vec, fc.vec, pv.vec, qv.vec)
      rownames(res) <- rownames(prop)
      
      if (mt.method == 'fdr') {
        res.final <- rbind(res.final, res[res[, 'Qvalue'] <= cutoff, , drop=F])
      }
      if (mt.method == 'raw') {
        res.final <- rbind(res.final, res[res[, 'Pvalue'] <= cutoff, , drop=F])
      }
    }
    
    if (!is.null(res.final)) {
      colnames(res.final) <- colnames(res)
    }
    return(list(pv.list=pv.list, fc.list=fc.list, qv.list=qv.list, m.list=m.list, res.final=res.final))
  }
  
}

visualize_differential_analysis <- function (data.obj, diff.obj,  grp.name=NULL, strata=NULL, test='Nonpara', mt.method='fdr', scale='sqrt', cutoff=0.15,
                                             taxa.levels=c('Phylum', 'Family', 'Genus'), ord=TRUE, eff.type='logP', indivplot=TRUE, colFnsC=NULL, colFnsF=NULL, subject=NULL,
                                             xsize=10, ann='', hei1=NULL, wid1=NULL, hei2=NULL, wid2=NULL) {
  
  # uniquefy names
  # For backward compatibility. Newer version will not need this and below. The old version has 'unclassified' which leads to duplicate names.
  # Newer version has 'Unclassified'. Case difference.
  
  # Check whether there is name duplication
  check.names <- NULL
  results <- list()
  obj0 <- diff.obj[[1]]
  for (level in names(obj0)) {
    obj <- obj0[[level]]
    # rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
    check.names <- c(check.names, rownames(obj))
  }
  
  if (sum(table(check.names) >= 2)) {
    data.obj <- uniquefy_taxa_names(data.obj)
    
    for (name1 in names(diff.obj)) {
      obj0 <- diff.obj[[name1]]
      for (level in names(obj0)) {
        obj <- obj0[[level]]
        # rownames(obj) <- gsub('unclassified', paste0('Unclassified',substr(level, 1, 1)), rownames(obj))
        rownames(obj) <- paste0(rownames(obj), substr(level, 1, 1))
        obj0[[level]] <- obj
      }
      diff.obj[[name1]] <- obj0
    }
    
  }
  
  fc.list <- diff.obj$fc.list
  qv.list <- diff.obj$qv.list
  pv.list <- diff.obj$pv.list
  if (test == 'Para') {
    fc.lc.list <- diff.obj$fc.lc.list
    fc.uc.list <- diff.obj$fc.uc.list
  }
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  
  ind <- !is.na(grp)
  data.obj <- subset_data(data.obj, ind)
  grp <- grp[ind]
  df <- df[ind, ]
  
  prop <- NULL
  eff <- eff.lc <- eff.uc <- NULL
  taxa.names <- NULL
  if (is.null(taxa.levels)) {
    LOIs <- names(qv.list)
  } else {
    LOIs <- taxa.levels
    if (sum(!(taxa.levels %in% names(qv.list)))) {
      stop('Taxa levels are not contained in differential abundance analysis results!\n')
    }
  }
  for (LOI in LOIs) {
    pv.vec <- pv.list[[LOI]]
    fc.vec <- fc.list[[LOI]]
    #qv.vec <- qvalue(pv.vec[, 1])$qvalues
    qv.vec <- qv.list[[LOI]]
    
    if (test == 'Para') {
      fc.lc.vec <- fc.lc.list[[LOI]]
      fc.uc.vec <- fc.uc.list[[LOI]]
    }
    
    if (mt.method == 'fdr') {
      taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
      taxa.name <- taxa.name[!is.na(taxa.name)]
    }
    
    if (mt.method == 'raw') {
      taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
      taxa.name <- taxa.name[!is.na(taxa.name)]
    }
    
    if (length(taxa.name) != 0) {
      prop0 <- data.obj$abund.list[[LOI]]
      prop0 <- t(t(prop0) / colSums(prop0))
      prop0 <-  prop0[taxa.name, , drop=F]
      if (ord == TRUE) {
        prop0 <- prop0[rev(order(rowMeans(prop0))), , drop=F]
      }
      prop <- rbind(prop, prop0)
      # currently using fold change
      if (test == 'Para') {
        eff <- rbind(eff, fc.vec[taxa.name, , drop=F])
        eff.lc <- rbind(eff.lc, fc.lc.vec[taxa.name, , drop=F])
        eff.uc <- rbind(eff.uc, fc.uc.vec[taxa.name, , drop=F])
      } else {
        if (eff.type == 'LFC') {
          eff <- c(eff, fc.vec[taxa.name, ])
        }
        if (eff.type == 'logP') {
          eff <- c(eff, sign(fc.vec[taxa.name, ]) * (-log10(pv.vec[taxa.name, ])))
        }
      }
      taxa.names <- c(taxa.names, taxa.name)
    }
  }
  results$taxa.names <- taxa.names
  if (length(taxa.names) == 0) {
    cat('No differential taxa! \n')
  } else {
    if (length(taxa.names) >= 2) {
      if (is.null(wid1) | is.null(hei1)) {
        wid1 <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
        hei1 <- 7
      }
      
      if (!is.null(grp.name)) {
        results$barplot_aggregate <- taxa_barplot_aggregate(prop, df, grp.name, strata, scale, xsize)
        results$boxplot_aggregate <- taxa_boxplot_aggregate(prop, df, grp.name, strata, scale, xsize) 
      }
      
      # currently fold change
      if (test == 'Para') {
        rownames(eff) <- rownames(eff.lc) <- rownames(eff.uc) <- taxa.names
        for (k in 1:ncol(eff)) {
          fold.dat.plot1 <- data.frame(Estimate=eff[, k], LCI=eff.lc[, k], UCI=eff.uc[, k], IV=taxa.names)
          results$effect_size <- plot_effect_size2(fold.dat.plot1)
        }
      } else {
        if (!is.na(eff[1])) {
          names(eff) <- taxa.names
          eff <- eff[!is.na(eff) & is.finite(eff)]
          eff <- sort(eff)
          taxa.names2 <- names(eff)
          if (is.null(wid2) | is.null(hei2)) {
            hei2 <- 4 + length(taxa.names2) / 20 * 3 
            wid2 <- 6
          }
          if (eff.type == 'LFC') {
            results$effect_size <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Log2 fold change')
          }
          if (eff.type == 'Spearman') {
            results$effect_size <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='Spearman correlation')
          }
          if (eff.type == 'logP') {
            results$effect_size <- plot_effect_size(taxa.names2, eff, levels(grp)[1], levels(grp)[2], ylab='-log10(P)')
          }					
        }
      }
      
      # create heatmp
      #			taxa.names2 <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
      if (!is.null(grp.name)) {
        
         results$prop_heatmap <- generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata),  ann=paste0(mt.method, '_', cutoff, '_', ann), colFnsC=colFnsC, colFnsF=colFnsF)
         results$rank_heatmap <-generate_taxa_heatmap(data.obj, taxa.levels='All', taxa=taxa.names, meta.info=c(grp.name, strata), data.type='R', ann=paste0(mt.method, '_', cutoff, '_Rank_', ann), colFnsC=colFnsC, colFnsF=colFnsF)
        try(
          results$biplot <- generate_taxa_biplot(data.obj, taxa=taxa.names, trans='sqrt', grp.name, ann=paste0(mt.method, '_', cutoff, '_', ann), varname.size = 1.5)	
        )	
      }
    }
    if (!is.null(grp.name)) {
      # Individual plots
      if (indivplot == TRUE) {
        results$taxa_boxplot <- generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=paste0(mt.method, '_', cutoff, '_', ann))
        if (!is.null(subject)) {
          generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, subject=subject, ann=paste0(mt.method, '_', cutoff, '_', ann, '_Paired'))
        }
        results$taxa_boxplot_binary <- generate_taxa_boxplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, scale='binary', ann=paste0(mt.method, '_', cutoff, '_', ann))
        results$taxa_barplot <- generate_taxa_barplot(data.obj, grp.name=grp.name, taxa.levels='All', strata = strata, taxa.name=taxa.names, ann=paste0(mt.method, '_', cutoff, '_', ann))
      }
      
    }
  }
  if(length(intersect(taxa.levels, c('Phylum', 'Family', 'Genus')))){
    results$cladogram <- make_cladogram(data.obj=data.obj, diff.obj=diff.obj, grp.name=grp.name, mt.method=mt.method, cutoff=cutoff)
  }
  return(results)
}

# Rev: 2018_03_08 add bootstrap standard error (normal approximation)
taxa_barplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='sqrt',  xsize=10, ylab='Proportion', error='ci', cutoff=0.00001) {
  
  if (scale == 'log') {
    prop[prop <= cutoff] <- cutoff
  }
  
  grp <- factor(df[, grp.name])
  
  if (is.null(strata)) {
    df2 <- data.frame(Group=grp, t(prop))
    df2 <- melt(df2)
    colnames(df2) <- c('Group', 'Taxa', 'Value')
    
    # Could be revised
    temp1 <- aggregate(Value ~ Group + Taxa, df2, mean)
    
    if (error == 'se') {
      temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) {
        sd(sapply(1:50, function(i) mean(sample(x, repl=TRUE))))
      })
    }
    if (error == 'ci') {
      temp2 <- aggregate(Value ~ Group + Taxa, df2, function(x) {
        2 * sd(sapply(1:50, function (i) mean(sample(x, repl=TRUE))))
      })
    }
    
    df2 <- cbind(temp1, temp2[, 3])
    colnames(df2) <- c('Group', 'Taxa', 'Mean', 'SE')
    
    limits <- aes(ymax = Mean + SE, ymin = ifelse(Mean - SE > 0, Mean - SE, 0))
    dodge <- position_dodge(width=0.90)
    
    obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) + 
      geom_bar(position=dodge, stat="identity", alpha=0.75) + 
      geom_bar(position=dodge, stat="identity", alpha=0.75, colour="black", show.legend=FALSE, size=0.25) +
      geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
      labs(y=paste(ylab), x='') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
      theme(legend.position="top", legend.title=element_blank()) + coord_flip()
    
    if (scale == 'sqrt') {
      obj1 <- obj1 + scale_y_sqrt(
        breaks = trans_breaks("sqrt", function(x) x^2),
        labels = trans_format("sqrt", math_format(.x^2)))
    }
    
    # To be revised
    if (scale == 'log') {
      obj1 <- obj1 + scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)))
    }
    if (scale == 'boxcox') {
      obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
    }
    
  } else {
    grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
    df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
    df2 <- melt(df2)
    colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')
    
    # Could be revised
    temp1 <- aggregate(Value ~ Group + Strata + Taxa, df2, mean)
    
    if (error == 'se') {
      temp2 <- aggregate(Value ~ Group + Strata +Taxa, df2, function(x) {
        #sd(x) / sqrt(length(x))
        sd(sapply(1:50, function(i) mean(sample(x, repl=TRUE))))
      })
    }
    if (error == 'ci') {
      temp2 <- aggregate(Value ~ Group + Strata + Taxa, df2, function(x) {
        #sd(x) / sqrt(length(x))
        2 * sd(sapply(1:50, function (i) mean(sample(x, repl=TRUE))))
      })
    }
    
    df2 <- cbind(temp1, temp2[, 4])
    colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Mean', 'SE')
    
    limits <- aes(ymax = Mean + SE, ymin = ifelse(Mean - SE > 0, Mean - SE, 0))
    dodge <- position_dodge(width=0.90)
    
    obj1 <- ggplot(df2, aes(x=Taxa, y=Mean, fill=Group)) + 
      geom_bar(position=dodge, stat="identity", alpha=0.75) + 
      geom_bar(position=dodge, stat="identity", alpha=0.75, colour="black", show.legend=FALSE, size=0.25) +
      geom_errorbar(limits, position=dodge, size=0.25, width=0.25) +
      facet_wrap(~Strata, ncol=1) +
      labs(y=paste(ylab), x='') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
      theme(legend.position="top", legend.title=element_blank())
    if (scale == 'sqrt') {
      obj1 <- obj1 + scale_y_sqrt(
        breaks = trans_breaks("sqrt", function(x) x^2),
        labels = trans_format("sqrt", math_format(.x^2)))
    }
    if (scale == 'log') {
      obj1 <- obj1 + scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)))
    }
    if (scale == 'boxcox') {
      obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
    }
  }
  return(obj1)
}

# Rev: 2018_03_08 Add geom_point and jsize
taxa_boxplot_aggregate <- function (prop, df, grp.name, strata=NULL, scale='none', xsize=10, jsize=NULL, ylab='Proportion', cutoff=0.00001) {
  grp <- factor(df[, grp.name])
  if (scale == 'log') {
    prop[prop <= cutoff] <- cutoff
  }
  if (is.null(jsize)) {
    nT <- nrow(prop) 
    jsize <- 2 / (nT %/% 20 + 1)
  }
  
  if (is.null(strata)) {
    df2 <- data.frame(Group=grp, t(prop))
    df2 <- melt(df2)
    colnames(df2) <- c('Group', 'Taxa', 'Value')
    
    dodge <- position_dodge(width=0.88)
    
    obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) + 
      geom_boxplot(position=dodge, alpha = ifelse(jsize == 0, 0.75, 0.75), outlier.alpha=ifelse(jsize == 0, 0.5, 0),  lwd=0.25, fatten=1)
    if (jsize != 0) {
      obj1 <- obj1 + geom_point(position=position_jitterdodge(dodge.width=0.88), size=jsize, alpha=0.3)
    }
    obj1 <- obj1 +
      labs(y=paste(ylab), x='') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
      theme(legend.position="top", legend.title=element_blank())
    if (scale == 'sqrt') {
      obj1 <- obj1 + scale_y_sqrt(
        breaks = trans_breaks("sqrt", function(x) x^2),
        labels = trans_format("sqrt", math_format(.x^2)))
    }
    if (scale == 'log') {
      obj1 <- obj1 + scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)))
    }
    if (scale == 'boxcox') {
      obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
    }
    
  } else {
    grp2 <- factor(paste(strata, df[, strata]), levels=paste(strata, levels(df[, strata])))
    df2 <- data.frame(Group=grp, Strata=grp2, t(prop))
    df2 <- melt(df2)
    colnames(df2) <- c('Group', 'Strata', 'Taxa', 'Value')
    
    dodge <- position_dodge(width=0.88)
    
    obj1 <- ggplot(df2, aes(x=Taxa, y=Value, fill=Group)) + 
      geom_boxplot(position=dodge, alpha = ifelse(jsize == 0, 0.75, 0.75),  outlier.alpha=ifelse(jsize == 0, 0.5, 0),  lwd=0.25, fatten=1) 
    if (jsize != 0) {
      obj1 <- obj1 + geom_point(position=position_jitterdodge(dodge.width=0.88), size=jsize, alpha=0.2)
    }
    obj1 <- obj1 +
      facet_wrap(~Strata, ncol=1) +
      labs(y=paste(ylab), x='') +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, size=xsize)) +
      theme(legend.position="top", legend.title=element_blank())
    if (scale == 'sqrt') {
      obj1 <- obj1 + scale_y_sqrt(
        breaks = trans_breaks("sqrt", function(x) x^2),
        labels = trans_format("sqrt", math_format(.x^2)))
    }
    # To be revised
    if (scale == 'log') {
      obj1 <- obj1 + scale_y_log10(
        breaks = trans_breaks("log10", function(x) 10^x),
        labels = trans_format("log10", math_format(10^.x)))
    }
    if (scale == 'boxcox') {
      obj1 <- obj1 + scale_y_continuous(breaks=c(0, 0.0001, 0.001, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 1.0), trans=boxcox_trans(1/3))
    }
    
  }
  return(obj1)
  
}

plot_effect_size <- function (month,  value, pos.lab, neg.lab, ylab, hjust1=1.3, hjust2=-0.3, lab.size=3, xsize=10) {
  month <- factor(month, levels=month)
  dtm <- data.frame(month=month, value=value)
  dtm$colour <- factor(ifelse(dtm$value < 0, neg.lab, pos.lab), levels=c(pos.lab, neg.lab))
  dtm$hjust <- ifelse(dtm$value > 0, hjust1, hjust2)
  obj <- ggplot(dtm, aes(month, value, label = month, hjust = hjust)) + 
    geom_text(aes(y = 0, colour = colour), size=lab.size) + 
    geom_bar(stat = "identity", aes(fill = colour)) +
    theme(axis.text.x = element_text(size=xsize)) +
    ylim(c(-max(abs(value))*1.1, max(abs(value))*1.1)) +
    coord_flip() + 
    scale_x_discrete(breaks = NULL) +
    labs(x = "", y = ylab) +
    theme(legend.position="top", legend.title=element_blank())
  return(obj)
}


plot_effect_size2 <- function (fold.dat.plot1, ylabel='log(Fold change)', is.ln=TRUE, ord=TRUE) {
  
  if (is.ln) {
    fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] <- fold.dat.plot1[, c('Estimate', 'LCI', 'UCI')] * log2(exp(1))
  }
  Alphas <- seq(1, 99, 2) / 100

  Multiplier <- seq(1, 0.01, len=50)
  zzTransparency <<- 1/(length(Multiplier)/4)
  
  fold.dat.plot1 <- data.frame(cbind(fold.dat.plot1, Scalar=rep(Multiplier, each = nrow(fold.dat.plot1))))
  fold.dat.plot1$Emphasis <- by(1 - seq(0, 1, length = length(Multiplier) + 1)[-(length(Multiplier) + 1)],
                                as.character(round(Multiplier, 5)), mean)[as.character(round(fold.dat.plot1$Scalar, 5))]
  
  fold.dat.plot1$IV <- factor(fold.dat.plot1$IV, unique(fold.dat.plot1$IV))
  
  if (ord) {
    fold.dat.plot1 <- fold.dat.plot1[order(as.character(fold.dat.plot1$IV)), ]
  }
  
  OutputPlot <- ggplot2::qplot(data = fold.dat.plot1, x = IV, y = Estimate,
                               ymin = Estimate - (Estimate -LCI)*Scalar, ymax = Estimate + (UCI - Estimate)*Scalar,
                               ylab = NULL, xlab = NULL, alpha = I(zzTransparency), colour = I(gray(0)), geom = "blank")
  
  OutputPlot <- OutputPlot + geom_hline(yintercept = 0, lwd = I(7/12), colour = I(hsv(0/12, 7/12, 7/12)), alpha = I(5/12))
  OutputPlot <- OutputPlot + geom_linerange(data = fold.dat.plot1, aes(size = 1/Emphasis), alpha = I(zzTransparency), colour = I(gray(0)))
  OutputPlot <- OutputPlot + scale_size_continuous() + guides(size=FALSE)
  OutputPlot <- OutputPlot + coord_flip() + geom_point(aes(x = IV, y = Estimate), colour = I(gray(0))) + theme_bw() + ylab(ylabel)
  return(OutputPlot)
}

# Rev: 2017_02_19 Add hei and wid
# Rev: 2018_01_15 Add qt.outlier=0.097
# Rev: 2018_03_06 position_jitter(h = 0), coord_cartesian
generate_taxa_boxplot <- function (data.obj,  grp.name, strata=NULL, scale='P',  taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'), 
                                   taxa.name='All', pseudo.ct=0.5, rm.outlier=T, qt.outlier=0.97, prev=0.1, minp=0.002, ann='All', subject=NULL, l.size=0.5, p.size=2.5, hei0=NULL, wid0=NULL) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  if (!is.null(subject)) {
    ID <- df[, subject]
  }
  obj <- NULL
  for (LOI in taxa.levels) {
    if (LOI == 'All') {
      if (taxa.name == 'All') stop('Taxa names are not required to be all when taxa level is also all!\n')
      headnames <- NULL
      prop <- NULL
      for (LOI2 in names(data.obj$abund.list)) {
        prop2 <- data.obj$abund.list[[LOI2]]
        taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
        if (length(taxa.name2) != 0) {
          if (scale == 'logP') {
            prop2 <- prop2 + pseudo.ct
          } 
          prop2 <- t(t(prop2) / colSums(prop2))			
          headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
          names(headnames2) <- rownames(prop2)
          prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
          headnames <- c(headnames, headnames2)
        }
      }
      
    } else {
      cat(LOI, "\n")
      prop <- data.obj$abund.list[[LOI]]
      if (scale == 'logP') {
        prop <- prop + pseudo.ct
      } 
      prop <- t(t(prop) / colSums(prop))
      
      headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
      names(headnames) <- rownames(prop)
      
      if (taxa.name == 'All') {
        prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
      } else {
        prop <- prop[taxa.name, , drop=FALSE]
      }
      
    }
    
    if (scale == 'logP') {
      prop <- log10(prop)
    }
    
    if (scale == 'sqrtP') {
      prop <- sqrt(prop)
    }
    
    if (scale == 'binary') {
      temp <- prop != 0
      prop[temp] <- 'Presence'
      prop[!temp] <- 'Absence'
    }
    if (is.null(hei0)) {
      hei <- 5
    } else {
      hei <- hei0
    }
    
    if (is.null(strata)) {
      if (is.null(wid0)) {
        wid <- 5
      } else {
        wid <- wid0
      }
      
    } else {
      if (is.null(wid0)) {
        wid <- 5.5
      } else {
        wid <- wid0
      }
    }
    
    if (scale == 'P') {
      ylab <- 'Proportion'
    } 
    if (scale == 'logP') {
      ylab <- 'log10(Proportion)'
    }
    if (scale == 'sqrtP') {
      ylab <- 'sqrt(Proportion)'
    }
    if (scale == 'binary') {
      ylab <- 'Count'
    } 
    
    if (is.null(strata)) {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        if (scale != 'binary') {
          if (scale == 'P') {
            if (rm.outlier == T) {
              ylims <- c(min(taxon.abund), quantile(taxon.abund, qt.outlier) * 1.25)
            } else {
              ylims <- range(taxon.abund)
            }
            
          } else {
            ylims <- range(taxon.abund)
          }
          
          if (is.null(subject)) {
            df2 <- data.frame(Value=taxon.abund, Group=grp)	
            dodge <- position_dodge(width=0.9)
            obj <- ggplot(df2, aes(x=Group, y=Value, col=Group)) +
              geom_boxplot(position=dodge, outlier.colour = NA, alpha=0.75) + 
              geom_jitter(alpha=0.6, size=3.0,  position = position_jitter(w = 0.1, h = 0)) +
              labs(y=ylab, title=headnames[taxon]) 
            if (rm.outlier) {
              obj <- obj + coord_cartesian(ylim = ylims)
            }
            #		ylim(ylims[1], ylims[2]) +
            obj <- obj + theme(legend.position="none")
          } else {						
            df2 <- data.frame(Value=taxon.abund, Group=grp, subject=ID)	
            dodge <- position_dodge(width=0.9)
            obj <- ggplot(df2, aes(x=Group, y=Value, shape=Group, group=subject)) +
              geom_point(size=p.size) +
              geom_line(size=l.size) +
              labs(y=ylab, title=headnames[taxon]) +
              #	ylim(ylims) +
              theme(legend.position="none")
          }
          
        } else {
          df2 <- data.frame(Value=taxon.abund, Group=grp)	
          obj <- ggplot(df2, aes(x=Group, fill=Value)) +
            geom_bar(width=.5) +
            labs(y=ylab, title=headnames[taxon]) +
            theme(legend.title=element_blank())
        }
      }	
    } else {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        if (scale != 'binary') {
          if (scale == 'P') {
            if (rm.outlier == T) {
              ylims <- c(min(taxon.abund), quantile(taxon.abund, qt.outlier) * 1.25)
            } else {
              ylims <- range(taxon.abund)
            }
            
          } else {
            ylims <- range(taxon.abund)
          }
          grp2 <- df[, strata]
          df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
          dodge <- position_dodge(width=0.9)
          obj <- ggplot(df2, aes(x=Strata, y=Value, col=Group, fill=Group)) +
            geom_boxplot(fill='white', alpha=0.6, position=dodge, outlier.colour = NA, alpha=0.75) + 
            geom_point(position=position_jitterdodge(dodge.width=0.9), size=3.0, alpha=0.6) +
            labs(y=ylab, x=strata, title=headnames[taxon])
          if (rm.outlier) {
            obj <- obj + coord_cartesian(ylim = ylims)
          }
          obj <- obj + theme(legend.position="none")
        } else {
          grp2 <- df[, strata]
          df2 <- data.frame(Value=taxon.abund, Group=grp, Strata=grp2)
          obj <- ggplot(df2, aes(x=Group, fill=Value)) +
            geom_bar(width=.5) +
            labs(y=ylab, title=headnames[taxon]) +
            facet_wrap(~ Strata) + 
            theme(legend.title=element_blank())
        }
      }
      
    }
    return(obj)
  }
  
}

# Add combined barplot with error bar and presence/absence bar (currently presence/absence bar is in generate_taxa_boxplot
# Rev: 2018_01_15 Add qt.outlier=0.97
# Rev: 2018_03_06 coord_cartesian
generate_taxa_barplot <- function (data.obj,  grp.name, strata=NULL, scale='P', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
                                   taxa.name='All', rm.outlier=T, qt.outlier=0.97, prev=0.1, minp=0.002, ann='All') {
  # To be completed
  df <- data.obj$meta.dat
  grp <- factor(df[, grp.name])
  obj <- NULL
  for (LOI in taxa.levels) {
    if (LOI == 'All') {
      if (taxa.name == 'All') stop('Taxa names are not required to be all when taxa level is also all!\n')
      headnames <- NULL
      prop <- NULL
      for (LOI2 in names(data.obj$abund.list)) {
        prop2 <- data.obj$abund.list[[LOI2]]
        taxa.name2 <- taxa.name[taxa.name %in% rownames(prop2)]
        if (length(taxa.name2) != 0) {
          if (scale == 'logP') {
            prop2 <- prop2 + 0.5
          } 
          prop2 <- t(t(prop2) / colSums(prop2))			
          headnames2 <- sapply(strsplit(rownames(prop2), ";"), paste, collapse="\n")
          names(headnames2) <- rownames(prop2)
          prop <- rbind(prop, prop2[taxa.name2, , drop=FALSE])
          headnames <- c(headnames, headnames2)
        }
      }
      
    } else {
      cat(LOI, "\n")
      prop <- data.obj$abund.list[[LOI]]
      if (scale == 'logP') {
        prop <- prop + 0.5
      } 
      prop <- t(t(prop) / colSums(prop))
      
      headnames <- sapply(strsplit(rownames(prop), ";"), paste, collapse="\n")
      names(headnames) <- rownames(prop)
      
      if (taxa.name == 'All') {
        prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
      } else {
        prop <- prop[taxa.name, , drop=FALSE]
      }
      
    }
    
    if (scale == 'logP') {
      prop <- log10(prop)
    }
    
    
    hei <- 5
    if (is.null(strata)) {
      wid <- 5
    } else {
      wid <- 4 * nlevels(df[, strata])
    }
    
    if (scale == 'P') {
      ylab <- 'Proportion'
    } else {
      ylab <- 'log10(Proportion)'
    }
    
    if (is.null(strata)) {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        taxon.abund2 <- taxon.abund
        if (scale == 'P') {
          if (rm.outlier == T) {
            ylims <- c(0, quantile(taxon.abund, qt.outlier) * 1.25)
            taxon.abund[taxon.abund > quantile(taxon.abund, qt.outlier) * 1.25] <- quantile(taxon.abund, qt.outlier) * 1.25
          } else {
            ylims <- c(0, max(taxon.abund))
          }
          
        } else {
          ylims <- range(taxon.abund)
        }
        df2 <- data.frame(x=factor(1:length(taxon.abund)), Value=taxon.abund, Group=grp)			
        dodge <- position_dodge(width=0.99)
        obj <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
          geom_bar(position=dodge, stat='identity', alpha=0.75) + 
          facet_grid(. ~ Group, scales='free_x', space="free") +
          labs(y=ylab, title=headnames[taxon])
        
        if (rm.outlier) {
          obj <- obj + coord_cartesian(ylim = ylims)
        }

        obj <- obj  + 
          xlab('') +
          theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
          theme(legend.position="none")
        mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, mean))
        obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint))		
        mean_df <- data.frame(Group=levels(grp), yint = tapply(taxon.abund2, grp, median))
        obj <- obj + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)

      }	
    } else {
      for (taxon in rownames(prop)) {
        taxon.abund <- prop[taxon, ]
        taxon.abund2 <- taxon.abund
        if (scale == 'P') {
          if (rm.outlier == T) {
            ylims <- c(0, quantile(taxon.abund, qt.outlier) * 1.25)
            taxon.abund[taxon.abund > quantile(taxon.abund, qt.outlier) * 1.25] <- quantile(taxon.abund, qt.outlier) * 1.25
          } else {
            ylims <- c(0, max(taxon.abund))
          }
          
        } else {
          ylims <- range(taxon.abund)
        }
        
        grp2 <- df[, strata]
        obj.list <- list()
        for (level in levels(grp2)) {
          ind <- grp2 %in% level
          df2 <- data.frame(x=factor(1:length(taxon.abund[ind])), Value=taxon.abund[ind], Group=grp[ind])	
          
          dodge <- position_dodge(width=0.99)
          obj0 <- ggplot(df2, aes(x=x, y=Value, fill=Group)) +
            geom_bar(position=dodge, stat='identity', alpha=0.75) + 
            facet_grid(. ~ Group, scales='free_x', space="free") +
            labs(y=ylab, title=headnames[taxon])
          
          if (rm.outlier) {
            obj0 <- obj0 + coord_cartesian(ylim = ylims)
          }
          
          
          obj0 <- obj0  + 
            xlab(paste(strata, level)) +
            theme(axis.ticks.x=element_blank(), axis.text.x = element_blank()) +
            theme(legend.position="none")
          mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], mean))
          obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint))		
          mean_df <- data.frame(Group=levels(grp[ind]), yint = tapply(taxon.abund2[ind], grp[ind], median))
          obj0 <- obj0 + geom_hline(data = mean_df, aes(yintercept = yint), linetype=2)
          obj.list[[level]] <- obj0		
        }
        obj <- multiplot(plotlist=obj.list, cols=nlevels(grp2))
      }	
      
    }
    return(obj)
  }		
}


generate_taxa_biplot <- function (data.obj, taxa, trans='sqrt', grp.name, ann='', ...) {
  
  grp <- data.obj$meta.dat[, grp.name]
  
  prop <- NULL
  for (LOI2 in names(data.obj$abund.list)) {
    ct <- data.obj$abund.list[[LOI2]]
    ct <- ct + 1
    prop0 <- t(t(ct) / colSums(ct))
    prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])	
  }
  colnames(prop) <- colnames(prop0)
  if (nrow(prop) != length(taxa)) {
    warnings('Some taxa not found in abundance lists! Please check the names!\n')
  }
  
  if (trans == 'normal') 	prop <- t(apply(prop, 1, function(x) qqnorm(x, plot=F)$x))
  if (trans == 'log') prop <- log(prop)
  if (trans == 'sqrt') prop <- sqrt(prop)
  if (trans == 'rank') prop <- t(apply(prop, 1, rank))
  
  wine.pca <- prcomp(t(prop), scale. = TRUE)
  g <- ggbiplot(wine.pca, obs.scale = 1, var.scale = 1, 
                groups = grp, ellipse = TRUE, circle = FALSE, ...) 
  g <- g + scale_color_discrete(name = '') + theme_bw()
  g <- g + theme(legend.direction = 'horizontal', legend.position = 'top') 
  return(g)
}

# Rev: 2017_02_19 key alignment modified, rowsep, colsep bug
# Rev: 2017_12_19 row colors fixed
generate_taxa_heatmap <- function (data.obj, taxa.levels='Genus', taxa='All', meta.info, sam.ord=NULL, data.type='P',  prev=0.1, minp=0.002, 
                                   row.col.dat='Phyla', phy.no=4, sepwidth=0.01, colsep=NULL, rowsep=NULL, sepcolor='black',
                                   white='white', colFnsC=NULL, colFnsF=NULL, Rowv=T, Colv=T, dendrogram='both', margins=c(5, 15), in.grid=F,  is.labCol=T, cexCol=1, cexRow=NULL,
                                   omas=c(1, 1, 1, 8), width=12, height=6, ann='All', return.obj=FALSE, ...) {
  df <- data.obj$meta.dat
  prop <- NULL
  if ('Species' %in% taxa.levels & !('Species' %in% names(data.obj$abund.list))) {
    data.obj$abund.list[['Species']] <- data.obj$otu.tab
    rownames(data.obj$abund.list[['Species']]) <- paste0("OTU", rownames(data.obj$otu.tab), ":", 
                                                         data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
  }
  
  for (LOI in taxa.levels) {
    cat(LOI, "\n")
    
    if (LOI == 'All') {
      if (taxa[1] == 'All') {
        stop("Please specify the taxa names that will be included in the heatmap!\n")
      } 
      prop <- NULL
      for (LOI2 in names(data.obj$abund.list)) {
        ct <- data.obj$abund.list[[LOI2]]
        prop0 <- t(t(ct) / colSums(ct))
        prop <- rbind(prop, prop0[intersect(rownames(prop0), taxa), , drop=FALSE])	
        
      }
      colnames(prop) <- colnames(prop0)
      if (nrow(prop) != length(taxa)) {
        warnings('Some taxa not found in abundance lists! Please check the names!\n')
      }
      
    } else {
      ct <- data.obj$abund.list[[LOI]]
      prop <- t(t(ct) / colSums(ct))
      
      if (taxa[1] != 'All') {
        prop <- prop[taxa, , drop=FALSE]
      } else {
        prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]	
      }			
    }
    
    # Deal with zeros
    if (data.type == 'B') {
      col.scheme <- c("lightyellow", "red")
      prop[, ] <- as.numeric(prop != 0)
      breaks <- c(-0.01, 0.01, 1.1)
    } 
    if (data.type == 'P'){
      col.scheme = c(white, brewer.pal(11, "Spectral"))
      ind.temp <- prop != 0
      minp <- min(prop[prop!=0])/1.1
      prop[prop==0] <- minp
      prop <- log10(prop)
      breaks <- c(log10(minp)-0.01, seq(log10(minp)+0.01, 0, len=12))
    }
    if (data.type == 'R'){
      col.scheme <- c(white, colorRampPalette(rev(brewer.pal(9, "RdYlBu")))(ncol(prop)-1))
      
      prop <- t(apply(prop, 1, function(x) {
        temp <- rank(x[x!=0])
        s <- (ncol(prop) - 1) / (max(temp) - min(temp))
        temp <- 1 + (temp - min(temp)) * s
        x[x!=0] <- temp
        x
      }))
      breaks <- seq(0, ncol(prop), len=ncol(prop)+1)
    }
  }
  
  phy <- sapply(strsplit(rownames(prop), ";"), function(x) x[1])
  
  obj <- main_heatmap(prop, colors=col.scheme) %>% 
    add_row_annotation(phy) %>%
    add_col_annotation(as.data.frame(data.obj$meta.dat[,meta.info, drop=FALSE])) %>% 
    add_row_clustering() %>% 
    add_col_clustering()
    
  return(obj)
  
}

clade.anno2 <- function(gtree, anno.data, alpha=0.2, anno.depth=3, anno.x=10, anno.y=40){
  short.labs <- c(letters, paste0(letters, letters))
  get_offset <- function(x) {(x*0.2+0.2)^2}
  get_angle <- function(node){
    data <- gtree$data
    sp <- ggtree:::get.offspring.df(data, node)
    sp2 <- c(sp, node)
    sp.df <- data[match(sp2, data$node),]
    mean(range(sp.df$angle))
  }
  anno.data = arrange(anno.data, node)
  hilight.color <- anno.data$color
  node_list <- anno.data$node
  node_ids <- (gtree$data %>% filter(label %in% node_list ) %>% arrange(label))$node
  anno <- rep('white', nrow(gtree$data))

  for(i in 1:length(node_ids)){
    n <- node_ids[i]
    color <- hilight.color[i]
    anno[n] <- color
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    offset <- get_offset(nodeClass)
    gtree <-
      gtree + geom_hilight(node=n, fill=color, alpha=alpha,
                           extend=offset)
  }
  gtree$layers <- rev(gtree$layers)
  gtree <- gtree + geom_point2(aes(size=I(nodeSize)), fill=anno, shape=21)
  ## add labels
  short.labs.anno <- NULL
  for(i in 1:length(node_ids)){
    n <- node_ids[i]
    mapping <- gtree$data %>% filter(node == n)
    nodeClass <- as.numeric(mapping$nodeClass)
    lab <- short.labs[1]
    short.labs <- short.labs[-1]
    color = anno[n]
    if(is.null(short.labs.anno)){
      short.labs.anno = data.frame(lab=lab, annot = mapping$label, node=n, color=color, stringsAsFactors = F)
    }else{
      short.labs.anno = rbind(short.labs.anno,
                              c(lab, mapping$label, n, color))
    }
    offset <- get_offset(nodeClass) - 0.4
    angle <- get_angle(n) + 90
    gtree <-
      gtree + geom_cladelabel(node=n, label=lab, fontsize=1.5+sqrt(nodeClass),
                              offset.text=offset, barsize=0, hjust=0.5)
  }
  anno_shapes = sapply(short.labs.anno$lab, utf8ToInt)
  labels = paste0(short.labs.anno$lab, ": ", gtree$data$label[node_ids])
  df2 <- unique(anno.data[,c("color", "cat")])
  df2$x=0
  df2$y=1
  df2$alpha=1
  df3 <- data.frame(color="white", cat="clade", x=0, y=0, alpha=1)
  df1 <- data.frame(color=anno[node_ids], cat=labels, x=0, y=1, alpha=alpha)
  df <- rbind(df2, df3, df1)
  p <- gtree + geom_rect(aes_(fill=~cat, xmin=~x, xmax=~x, ymin=~y, ymax=~y), data=df, inherit.aes=F) + guides(fill=guide_legend(override.aes=list(fill=alpha(df$color, df$alpha)))) + theme(legend.position="right")
  p
}

make_cladogram <- function(data.obj, diff.obj, grp.name, cutoff=0.05, prev=0.1, minp=0.002, mt.method="fdr"){
  
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  levels(grp) <- paste0(1:nlevels(grp), levels(grp))

  qv.list <- diff.obj$qv.list
  pv.list <- diff.obj$pv.list
  fc.list <- diff.obj$fc.list
  
  otu.name.12 <- data.obj$otu.name
  otu.tab.12 <- data.obj$otu.tab
  tax.family.a <- NULL
  for (i in 1:6) {
    tax.family <- apply(otu.name.12[, 1:i, drop=F], 1, paste, collapse=".")
    phlan.tab <- aggregate(otu.tab.12, by=list(tax.family), FUN=sum)
    rownames(phlan.tab) <- phlan.tab[, 1]
    phlan.tab <- as.matrix(phlan.tab [, -1, drop=F])
    phlan.tab <- t(t(phlan.tab) / colSums(phlan.tab))	
    phlan.tab <- phlan.tab[rowSums(phlan.tab!=0) > prev*ncol(phlan.tab) & rowMaxs(phlan.tab) > minp, , drop=F]
    tax.family.a <- c(tax.family.a, rownames(phlan.tab))
  }
  alias.a <- sapply(strsplit(tax.family.a, "\\."), function(x) {
    if (length(x) >= 2) {
      if (length(x) == 2) {
        return(x[2])
      } else {
        return(paste0(x[2], ";", x[length(x)]))
      }
    } else {
      return(x)
    }
  })
  
  taxa.names <- NULL
  abundant.grp.names <- NULL
  
  for (LOI in setdiff(names(qv.list), 'Species')) {
    qv.vec <- qv.list[[LOI]]
    pv.vec <- pv.list[[LOI]]
    
    fc.vec <- fc.list[[LOI]]
    if (mt.method == 'fdr') {
      taxa.name <- rownames(qv.vec)[qv.vec <= cutoff]
      abundant.grp.name <- apply(fc.vec, 1, function(x) levels(grp)[as.numeric(x >= 0) + 1])[qv.vec <= cutoff]
    }
    
    if (mt.method == 'raw') {
      taxa.name <- rownames(pv.vec)[pv.vec <= cutoff]
      abundant.grp.name <- apply(fc.vec, 1, function(x) levels(grp)[as.numeric(x >= 0) + 1])[pv.vec <= cutoff]
    }
    
    if (length(taxa.name) != 0) {
      taxa.names <- c(taxa.names, taxa.name)
      abundant.grp.names <- c(abundant.grp.names, abundant.grp.name)
    }
  }
  
  # remove 'unclassified'
  abundant.grp.names <- abundant.grp.names[!grepl('unclassified', taxa.names, ignore.case=T)]
  taxa.names <- taxa.names[!grepl('unclassified', taxa.names, ignore.case=T)]
  
  ind <- match(taxa.names, alias.a)
  na.ind <- !is.na(ind)
  
  phlan <- cbind(tax.family.a, '1.5', '\t', '\t', '-')
  
  phlan[ind[na.ind], 2:5] <- cbind('3.5',  abundant.grp.names[na.ind], '3.0', '-')
  
  tax_chars2 <- c('k', 'p', 'c', 'o', 'f', 'g')
  
  phlan2 <- phlan
  
  phlan2[,1] <- unlist(lapply(phlan[,1], function(x){
    spl <- strsplit(x, '\\.')
    paste(paste0(tax_chars2, "__", spl[[1]])[1:length(spl[[1]])], collapse='|')
  }))
  
  tips <- names(tax.family[which(tax.family %in% tax.family.a)])
  tree <- data.obj$tree
  pruned <- drop.tip(data.obj$tree,data.obj$tree$tip.label[-match(tips,data.obj$tree$tip.label)])
  pruned$tip.label <- tax.family[pruned$tip.label]
  
  nodes <- lapply(unique(pruned$tip.label[duplicated(pruned$tip.label)]), function(x){
    ggtree::MRCA(pruned, tip=c(x,x))
  })
  
  pruned2 <- ggtree(pruned, layout="circular")
  
  for(x in unlist(nodes)){
    pruned2 <- ggtree::collapse(pruned2, node=x)
  }
  
  ###Nodes to highlight
  highlight <- sapply(strsplit(phlan2[which(phlan2[,2] == "3.5"),][,1], "\\|"), tail, 1)
  nolight <- sapply(strsplit(phlan2[which(phlan2[,2] != "3.5"),][,1], "\\|"), tail, 1)
  cat <- phlan2[which(phlan2[,2] == "3.5"),3]
  nocat <- phlan2[which(phlan2[,2] != "3.5"),3]
  nocat[nocat=="\t"] <- "NA"
  anno.data <- data.frame(node=highlight, cat=cat)
  palette <- I(brewer.pal(max(nlevels(anno.data$cat), 3), name='Set1'))[1:nlevels(anno.data$cat)]
  anno.data$color <- palette[anno.data$cat]
  testtree <- parseLefse(phlan)
  
  p <- tree.backbone(testtree)
  p <- clade.anno2(p, anno.data)
  p
}

parseLefse <- function(phlan, index=1, header=FALSE, delimiter='\\.', node.size.scale=1, node.size.offset=1){
  taxtab2 <- as.data.frame(phlan[,1:3])
  names(taxtab2)[1] <- 'tax'
  names(taxtab2)[2] <- 'score'
  names(taxtab2)[3] <- 'grp'
  taxtab2$tax <- as.character(taxtab2$tax)
  taxtab2 <- taxtab2 %>% slice(grep('unclassified', .[,index], invert=TRUE))
  
  tax_chars2 <- c('k', 'p', 'c', 'o', 'f', 'g')
  tax_split2 <- lapply(taxtab2$tax, function(x){
    spl <- strsplit(x, '\\.')
    tmp <- paste0(tax_chars2, "__", spl[[1]])[1:length(spl[[1]])]
    gsub("-","",tmp)
  })     ## split into different taxonomy levels
  child2 <- vapply(tax_split2, tail, n=1, '')
  tax_class2 <- do.call(rbind, strsplit(child2, '__'))[,1]
  parent2 <- vapply(tax_split2, function(x) ifelse(length(x)>1, x[length(x)-1], 'root'), '')
  isTip2 <- !child2 %in% parent2
  index2 <- c()
  index2[isTip2] <- 1:sum(isTip2)
  index2[!isTip2] <- (sum(isTip2)+1):length(isTip2)
  ## tips comes first
  mapping2 <- data.frame(node=index2, row.names=make.names(child2, unique=TRUE), taxa.names=child2, isTip2, check.names = FALSE)
  edges2 <- cbind(mapping2[parent2,]$node, mapping2$node)
  edges2 <- edges2[!is.na(edges2[,1]),]
  
  a <- node.size.scale
  b <- node.size.offset
  mapping2$nodeSize <- a+b
  ###mapping2$nodeSize <- a*log(mapping$taxaAbun) + b
  mapping2$nodeClass <- factor(tax_class2, levels = rev(tax_chars2))
  
  mapping2 <- mapping2[order(mapping2$node),]
  
  node.label2 <- rownames(mapping2)[!mapping2$isTip2]
  phylo <- structure(list(edge = edges2,
                          node.label = node.label2,
                          tip.label = rownames(mapping2)[mapping2$isTip2],
                          edge.length=rep(1, nrow(edges2)),
                          Nnode = length(node.label2)
  ),
  class = "phylo")
  
  d <- mapping2 %>% select_(~-isTip2)
  testtree <- treedata(phylo = phylo, data = as_data_frame(d))
}

predictionRF <- function (data.obj,  resp.name, formula=NULL, taxa.level='Species', binary=FALSE, prev=0.1, minp=0.002, B=50, seed=123, 
                          boruta.level=c('Confirmed', 'Tentative'), ann='',...) {

  response <- data.obj$meta.dat[, resp.name]
  taxa.names <- NULL
  if (!is.null(formula)) {
    adj <- model.matrix(as.formula(formula), data.obj$meta.dat)
    adj.var <- colnames(adj)
  } 
  
  if (taxa.level == 'Species') {
    if (taxa.level %in% names(data.obj$abund.list)) {
      ct <- data.obj$abund.list[[taxa.level]]
    } else {
      # Accomodate different version
      ct <- data.obj$otu.tab
      rownames(ct) <- paste0("OTU", rownames(ct), ":", data.obj$otu.name[, 'Phylum'], ";", data.obj$otu.name[, 'Genus'])
      data.obj$abund.list[['Species']] <- ct
    }
  } else {
    ct <- data.obj$abund.list[[taxa.level]]
  }
  
  prop <- t(t(ct) / colSums(ct))
  prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
  prop <- t(prop)
  
  if (binary == TRUE) {
    prop <- (prop != 0)
  }
  
  original.names <- colnames(prop)
  set.seed(seed)
  
  if (is.null(formula)) {
    performance <- matrix(0, nrow=B, ncol=2)
    colnames(performance) <- c("RF_M","Guess")
    roc.list <- list(RF_M=NULL)
    lab.list <- roc.list
  } else {
    performance <- matrix(0, nrow=B, ncol=4)
    colnames(performance) <- c("RF_M", "RF_CF", "RF_M+CF", "Guess")
    roc.list <- list("RF_M"=NULL, "RF_CF"=NULL, "RF_M+CF"=NULL)
    lab.list <- roc.list
  }
  
  
  colnames(prop) <- gsub(";", "_", colnames(prop))
  colnames(prop) <- gsub(":", "_", colnames(prop))
  colnames(prop) <- gsub("-", "_", colnames(prop))
  colnames(prop) <- gsub("\\.", "_", colnames(prop))
  names(original.names) <- colnames(prop)
  if (!is.null(formula)) {
    names(adj.var) <- adj.var
    original.names <- c(original.names, adj.var)
  }
  
  
  if (!is.null(formula)) {
    padj <- cbind(prop, adj)
  } else {
    padj <- prop
  }
  
  I <- nrow(prop)
  cat('Begin to bootstrap ...\n')
  for(b in 1:B){
    if (b %% 50 == 0) cat (".")
    err <- try({
      bsample <- sample(1:I, I, replace=T)
      if (is.factor(response)) {
        rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
        performance[b, "RF_M"] <- mean(response[-bsample]!= rf1$test$predicted)
        performance[b, "Guess"]<- mean(response[-bsample] != levels(response)[which.max(tabulate(response[bsample]))])
        roc.list[['RF_M']] <- c(roc.list[['RF_M']], rf1$test$vote[, 1])
        lab.list[['RF_M']] <- c(lab.list[['RF_M']], response[-bsample])
        if (!is.null(formula)) {
          rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
          rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample], sampsize=table(response[bsample]), ...)
          performance[b, "RF_CF"] <- mean(response[-bsample]!= rf2$test$predicted)
          performance[b, "RF_M+CF"] <- mean(response[-bsample]!= rf3$test$predicted)
          roc.list[['RF_CF']] <- c(roc.list[['RF_CF']], rf2$test$vote[, 1])
          roc.list[['RF_M+CF']] <- c(roc.list[['RF_M+CF']], rf3$test$vote[, 1])
          lab.list[['RF_CF']] <- c(lab.list[['RF_CF']], response[-bsample])
          lab.list[['RF_M+CF']] <- c(lab.list[['RF_M+CF']], response[-bsample])
        }
        
      } else {
        rf1 <- randomForest(x=prop[bsample, ], y=response[bsample], xtest=prop[-bsample, ], ytest=response[-bsample], ...)
        performance[b, "RF_M"] <- mean((response[-bsample] - rf1$test$predicted)^2)
        performance[b, "Guess"]<- mean((response[-bsample] - mean(response[bsample]))^2)
        
        if (!is.null(formula)) {
          
          rf2 <- randomForest(x=adj[bsample, ], y=response[bsample], xtest=adj[-bsample, ], ytest=response[-bsample],  ...)
          rf3 <- randomForest(x=padj[bsample, ], y=response[bsample], xtest=padj[-bsample, ], ytest=response[-bsample],  ...)
          performance[b, "RF_CF"] <- mean((response[-bsample] - rf2$test$predicted)^2)
          performance[b, "RF_M+CF"] <- mean((response[-bsample]- rf3$test$predicted)^2)
        }
      }
    })
    if (inherits(err, "try-error")) {
      next
    }
    
  }
  
  fridman <- NULL
  if (is.null(formula)) {
    fridman <- cat("Fridman.test p value: ", friedman.test(performance)$p.value, "\n")
  } else {
    fridman <- cat("Fridman.test p value (M+CF vs CF): ", friedman.test(performance[, c('RF_M+CF', 'RF_CF')])$p.value, "\n")
  }

  
  classification_error <- NULL
  if (!is.null(formula)) {
    performance2 <- performance[, c('Guess', 'RF_CF', 'RF_M', 'RF_M+CF')]
  } else {
    performance2 <- performance
  }
  if (is.factor(response)) {
    classification_error <- ggplot(melt(as.data.frame(performance2)), aes(x=variable, y=value)) + 
                              geom_boxplot(aes(fill=variable)) + 
                              guides(fill=FALSE) + 
                              theme(text = element_text(size=20), axis.title.x=element_blank()) + 
                              ylab("Classification error")
  } else {
    classification_error <-ggplot(melt(as.data.frame(performance2)), aes(x=variable, y=value)) + 
                              geom_boxplot(aes(fill=variable)) + 
                              guides(fill=FALSE) + 
                              theme(text = element_text(size=20), axis.title.x=element_blank()) +
                              ylab("PMSE")
  }	

  
  # ROC curve - only compare to the first level of the factor 
  if (is.factor(response)) {
    lab.list <- lapply(lab.list, function(x) {
      y <- factor(x)
      levels(y) <- levels(response)
      y
    }
    )
    roc <- createROC(roc.list, lab.list, pos.lab=levels(response)[1], file.name=NULL)
  }
  
  if (is.factor(response)) {
    rf <- randomForest(x=padj, y=response, sampsize=table(response), importance=T, ...)
  } else {
    rf <- randomForest(x=padj, y=response,  importance=T, ...)
  }
  
  importance <- as.data.frame(rf$importance)
  
  mean_dec_acc <- NULL
  mean_dec_gini <- NULL
  inc_mse <- NULL
  inc_node_purity <- NULL
  if (is.factor(response)) {
    mean_dec_acc <- importance[rev(order(importance$MeanDecreaseAccuracy)), ]
    mean_dec_gini <- importance[rev(order(importance$MeanDecreaseGini)), ]
  } else {
    inc_mse <- importance[rev(order(importance[, '%IncMSE'])), ]
    inc_node_purity <- importance[rev(order(importance$IncNodePurity)), ]
  }
  
  
  # Boruta Feature Selection - All subset selection (can't solve the confounding problem)
  obj.Boruta <- Boruta(padj, response, doTrace = 2, num.threads=15)	
  feat_select_tab <- obj.Boruta$finalDecision
  
  pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_", taxa.level, '_', ann,  ".pdf"), height=6, width=10)
  par(mar=par('mar') + c(3, 0, 0, 0))
  plot(obj.Boruta, main = "Feature selection by Boruta", ylab="Importance z-score", lwd = 0.5, las = 3, xlab = "",
       cex=1 / (ncol(prop)/50), cex.axis=0.25*200/ncol(prop), yaxt='n')
  axis(2, cex.axis=1)
  dev.off()
  #sink()
  
  pdf(paste0("Taxa_Random_forest_Boruta_Feature_Selection_Significant_Only_", taxa.level, '_', ann, ".pdf"), height=6, width=6)
  par(mar=par('mar') + c(3, 0, 0, 0))
  otu.ids <- plot.Boruta2(obj.Boruta, main = "Feature selection by Boruta", lwd = 0.5, las = 3, 
                          ylab="Importance z-score", xlab = "", cex.axis=1, yaxt='n')
  axis(2, cex.axis=1)
  dev.off()
  
  if ('Confirmed' %in% boruta.level) {
    taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed')]] 
    if (!is.null(formula)) {
      taxa.names <- setdiff(taxa.names, adj.var)
      aug.var <- adj.var
    } else {
      aug.var <- NULL
    }
    
    if (length(taxa.names) > 0) {
      if (is.factor(response)) {
        generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
        generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        if (length(taxa.names) > 1) {
          generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
          generate_taxa_boxplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
          build.decision.tree(data.obj,  resp.name=resp.name, aug.var=aug.var, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann)) 
        }
        gen <- data.obj$abund.list[[taxa.level]]
        gen <- t(t(gen) / colSums(gen))
        dat <- t(gen[taxa.names, ,drop=FALSE])
        colnames(dat) <- paste0('V', 1:ncol(dat))
        response2 <- factor(-as.numeric(response)+2)
        
        mylr <- function(formula, train, test){
          model <- randomForest(formula, data=train)
          x <- predict(model, newdata=test, type='prob')[, 2]
        }
        
        ACC <- Daim(response2 ~., model=mylr, data=as.data.frame(dat), labpos="1",
                    control=Daim.control(method="boot", number=100), cutoff="0.632+")
        pdf(paste0('BorutaFeatures_Confirmed_ROC_', taxa.level, '_0.632+', ann, '.pdf'), height=6, width=6)
        plot(ACC, main='Boruta taxa', method='0.632+', lwd=1.5)	
        abline(0, 1, col='black')
        x <- ACC
        legend("bottomright", legend = paste("AUC:", 
                                             formatC(Daim::auc(x$roc$sens632p, x$roc$spec632p), 
                                                     digits = max(3, getOption("digits") - 3))),
               inset = 0.01)
        dev.off()	
      } else {
        data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
        levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
        generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Confirmed_', ann))
        generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        if (length(taxa.names) > 1) {
          generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
          generate_taxa_boxplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Confirmed_', ann))
        }
      }
    }
  }
  
  if ('Tentative' %in% boruta.level) {
    taxa.names <- original.names[names(obj.Boruta$finalDecision)[obj.Boruta$finalDecision %in% c('Confirmed',  'Tentative')]] 
    if (!is.null(formula)) {
      taxa.names <- setdiff(taxa.names, adj.var)
      aug.var <- adj.var
    } else {
      aug.var <- NULL
    }
    
    if (length(taxa.names) > 0) {
      
      if (is.factor(response)) {		
        generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        generate_taxa_boxplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
        generate_taxa_barplot(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        if (length(taxa.names) > 1) {
          generate_taxa_barplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
          generate_taxa_boxplot_aggregate(data.obj, grp.name=resp.name, taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
          build.decision.tree(data.obj,  resp.name=resp.name, aug.var=aug.var, taxa.level=taxa.level, binary=binary, taxa=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann)) 
        }
        
        gen <- data.obj$abund.list[[taxa.level]]
        gen <- t(t(gen) / colSums(gen))
        dat <- t(gen[taxa.names, ,drop=FALSE])
        colnames(dat) <- paste0('V', 1:ncol(dat))
        response2 <- factor(-as.numeric(response)+2)
        
        mylr <- function(formula, train, test){
          model <- randomForest(formula, data=train)
          x <- predict(model, newdata=test, type='prob')[, 2]
        }
        
        ACC <- Daim(response2 ~., model=mylr, data=as.data.frame(dat), labpos="1",
                    control=Daim.control(method="boot", number=100), cutoff="0.632+")
        pdf(paste0('BorutaFeatures_Tentative_ROC_', taxa.level, '_0.632+', ann, '.pdf'), height=6, width=6)
        plot(ACC, main='Boruta taxa', method='0.632+', lwd=2)	
        abline(0, 1, col='black')
        x <- ACC
        legend("bottomright", legend = paste("AUC:", 
                                             formatC(Daim::auc(x$roc$sens632p, x$roc$spec632p), 
                                                     digits = max(3, getOption("digits") - 3))),
               inset = 0.01)
        dev.off()	
        
      } else {
        data.obj$meta.dat[, paste0(resp.name, '_b')] <- factor(data.obj$meta.dat[, resp.name] < median(data.obj$meta.dat[, resp.name]))
        levels(data.obj$meta.dat[, paste0(resp.name, '_b')]) <- c('High', 'Low')
        generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        generate_taxa_boxplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, scale='binary', ann=paste0('BorutaFeatures_Tentative_', ann))
        generate_taxa_barplot(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        if (length(taxa.names) > 1) {
          generate_taxa_barplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
          generate_taxa_boxplot_aggregate(data.obj, grp.name=paste0(resp.name, '_b'), taxa.levels=taxa.level, taxa.name=taxa.names, ann=paste0('BorutaFeatures_Tentative_', ann))
        }
      }
    }
  }
  return(list(fridman=fridman, classification_error = classification_error, roc=roc, mean_dec_acc = mean_dec_acc,
              mean_dec_gini = mean_dec_gini, inc_mse = inc_mse, inc_node_purity = inc_node_purity, feat_select_tab = feat_select_tab, taxa.names=taxa.names))
}

createROC <- function (pv.list, lab.list, pos.lab='1', file.name='ROC.pdf', width = 6, height = 6) {
  require(ROCR)
  n <- length(pv.list)
  aucs <- numeric(n)
  names(aucs) <- names(pv.list)
  roc <- NULL
  #	cols <- scales::hue_pal()(n)
  cols <- rep(c('red', 'blue', 'orange', 'cyan', 'purple'), ceiling(n/5))[1:n]
  ltys <- rep(c(1, 2), ceiling(n/2))[1:n]
  if(!is.null(file.name)){
    pdf(file.name, height=6, width=6)
  }
  for (i in 1:n) {
    
    cat("*")
    pv.mat <- pv.list[[i]]
    lab.mat <- lab.list[[i]]
    
    pred <- prediction(pv.mat, lab.mat==pos.lab)
    perf <- performance(pred, "tpr", "fpr")
    aucs[i] <- mean(unlist(performance(pred, 'auc')@y.values))
    df <- data.frame(x=unlist(perf@x.values), y=unlist(perf@y.values))
    roc <- ggplot(df, aes(x,y)) + 
      geom_line() + 
      geom_abline(intercept=0, slope=1, linetype="dashed") + 
      xlab("Average false positive rate") + 
      ylab("Average true positive rate") + 
      annotate("text", x=Inf, y=-Inf, hjust=1, vjust=0, label=paste0(names(pv.list), "(AUC:", round(aucs, 3), ")"), parse=TRUE, size=5)
  }
  
  if(!is.null(file.name)){
    dev.off()
  }
  return(roc)
}

plot.Boruta2 <- function (x, colCode = c("green", "yellow", "red", "blue"), sort = TRUE,
                          whichShadow = c(TRUE, TRUE, TRUE), col = NULL, xlab = "Attributes", ids=NULL,
                          ylab = "Importance", ...)
{
  if (class(x) != "Boruta")
    stop("This function needs Boruta object as an argument.")
  lz <- lapply(1:ncol(x$ImpHistory), function(i) x$ImpHistory[is.finite(x$ImpHistory[, i]), i])
  names(lz) <- colnames(x$ImpHistory)
  numShadow <- sum(whichShadow)
  lz <- lz[c(rep(TRUE, length(x$finalDecision)), whichShadow)]
  col <- Boruta:::generateCol(x, colCode, col, numShadow)
  if (sort) {
    ii <- order(sapply(lz, median))
    lz <- lz[ii]
    col <- col[ii]
  }
  names(lz) <- gsub("^X",  "", names(lz))
  if (is.null(ids)) {
    len <- sum(x$finalDecision %in% c('Confirmed', 'Tentative'))
    ind <- (length(lz) - len) : length(lz)
  } else {
    ind <- match(ids, names(lz))
  }
  
  boxplot(lz[ind], xlab = xlab, ylab = ylab, col = col[ind], ...)
  invisible(x)
  names(lz[ind])
}

# Rev: 2017_05_19: new formula for width
#  Ad presence and absence bar
generate_taxa_barplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
                                             taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8)) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- df[, grp.name]
  plotlist <- list()
  for (i in 1:length(taxa.levels)) {
    LOI <- taxa.levels[i]
    cat(LOI, "\n")
    prop <- data.obj$abund.list[[LOI]]
    prop <- t(t(prop) / colSums(prop))
    
    if (taxa.name == 'All') {
      prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
    } else {
      prop <- prop[taxa.name, , drop=FALSE]
    }
    
    # decreasing
    prop <- prop[rev(order(rowMeans(prop))), ]
    
    if (is.null(wids) | is.null(heis)) {
      wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
      wid <- sqrt(nlevels(grp) / 2) * wid
      hei <- 7
    } else {
      wid <- wids[i]
      hei <- heis[i]
    }
    
    
    plotlist[[i]] <- taxa_barplot_aggregate (prop, df, grp.name, strata, scale, xsize[i]) 
  }
  return(plotlist)
}

generate_taxa_boxplot_aggregate <- function (data.obj,  grp.name, strata=NULL, scale='sqrt', taxa.levels=c('Phylum', 'Class', 'Order', 'Family', 'Genus'),
                                             taxa.name='All',  prev=0.1, minp=0.002, ann='All', wids=NULL, heis=NULL, xsize=c(14, 12, 10, 9, 8, 6, 4), jsize = 0) {
  # To be completed
  df <- data.obj$meta.dat
  grp <- factor(df[, grp.name])
  plotlist <- list()
  for (i in 1:length(taxa.levels)) {
    LOI <- taxa.levels[i]
    cat(LOI, "\n")
    prop <- data.obj$abund.list[[LOI]]
    prop <- t(t(prop) / colSums(prop))
    
    if (taxa.name == 'All') {
      prop <- prop[rowMaxs(prop) > minp & rowSums(prop!=0) > prev*ncol(prop), , drop=FALSE]
    } else {
      prop <- prop[taxa.name, , drop=FALSE]
    }
    
    # decreasing
    prop <- prop[rev(order(rowMeans(prop))), ]
    
    if (is.null(wids) | is.null(heis)) {
      wid <- 7 * ifelse(nrow(prop) / 30 < 1, 1, nrow(prop) / 30)
      wid <- sqrt(nlevels(grp) / 2) * wid
      
      hei <- 7
    } else {
      wid <- wids[i]
      hei <- heis[i]
    }
    
    
    plotlist[[i]] <- taxa_boxplot_aggregate (prop, df, grp.name, strata, scale, xsize[i], jsize=jsize) 
  }
  return(plotlist)
}

build.decision.tree <- function(data.obj,  resp.name, taxa.level='Species', binary=FALSE, taxa, aug.var=NULL, ann='All') {
  ann <- paste(taxa.level, ann, sep="_")
  response <- data.obj$meta.dat[, resp.name]
  
  ct <- data.obj$abund.list[[taxa.level]]
  prop <- t(t(ct) / colSums(ct))
  prop <- prop[taxa, , drop=F] 	
  if (binary == TRUE) {
    prop <- (prop != 0)
  }
  
  dat <- as.data.frame(t(prop))
  # Rev: 2017_02_17 Add additional variables from meta dat
  if (!is.null(aug.var)) {
    dat <- cbind(dat, data.obj$meta.dat[, aug.var])
  }
  
  dat <- data.frame(dat, response)
  try(
    if (is.factor(response)) {
      fit <- rpart(response ~ ., method="class", data=dat)
      post(fit, file = paste0("Taxa_Unpruned_Classification_tree_", ann, ".ps"), title = "Unpruned Classification Tree")
      pfit<- prune(fit, cp= fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
      post(pfit, file = paste0("Taxa_Pruned_Classification_", ann, ".ps"), title = "Pruned Classification Tree")
      
    } else {
      fit <- rpart(response ~ ., method="anova", data=dat)
      post(fit, file = paste0("Taxa_Unpruned_Regression_tree_", ann, ".ps"), title = "Unpruned Regression Tree")
      pfit<- prune(fit, cp= fit$cptable[which.min(fit$cptable[,"xerror"]),"CP"])
      post(pfit, file = paste0("Taxa_Pruned_Regression_tree_", ann, ".ps"), title = "Pruned RegressionTree")
    }
  )	
}
