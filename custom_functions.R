# customized functions
filter_dmrs <- function(dmrs,
                        pval = 0.05,
                        fwer = 0.05,
                        min_abs_value = 0.5) {
  #####
  # Function for filtering the bumphunter results
  # dmrs = raw result of bumphunter
  # pval = maximum p-value (here: p-value-area)
  # fwer = maximum fwer (here: fwer-area)
  # min_abs_value = minimal value for the mean of the "value" (=average effect size)
  #####
  
  bh_filt <- dplyr::filter(
    dmrs,
    .data[["p.value"]] <= .env$pval,
    .data[["fwer"]]    <= .env$fwer,
    #abs(.data[["value"]])  >= .env$min_abs_value
  )
  bh_filt %>% 
    arrange(bh_filt$value)
}

makeHeatmap_dmr <- function(pdfname, dmrs, anno, myNorm, targets, pd_col, compare_values) {
  ##############
  # Function for creating heatmaps for DMRs
  # pdfname = Filename for pdf
  # dmrs = Dataframe of annotated dmrs (limited 
  # anno = Annotation for the array
  # myNorm = beta-values 
  # targets = Sample Sheet (attributes)
  # pd_col = string of column to be evalutated
  # compare_values = list of strings of values of the attribute to be evaluated
  ##############
  
  
  pdf(pdfname, width = 10, height = 8)
  
  for (i in 1:nrow(dmrs)) {
    #Aktuelles DMR
    dmr_row <- dmrs[i, ]
    dmr_chr <-  paste0("chr", gsub("^chr", "", dmr_row$seqnames))
    dmr_start <- dmr_row$start
    dmr_end <- dmr_row$end
    
    
    cat("Starte Heatmap für DMR", i, "\n")
    
    #CpGs im DMR
    dmr_probes <- anno$Name[
      anno$chr == dmr_chr &
        anno$pos >= dmr_start &
        anno$pos <= dmr_end
    ]
    
    common_probes <- intersect(dmr_probes, rownames(myNorm))
    
    
    common_gene_probes_all <- common_probes
    
    if(length(common_gene_probes_all) == 0) {
      cat("Kein CpG für DMR ", i, "gefunden – übersprungen.\n")
      next
    }
    
    #Beta-Werte extrahieren & nach Position sortieren
    gene_beta <- myNorm[common_gene_probes_all, , drop = FALSE]
    
    # CpGs nach genomischer Position sortieren
    probe_anno <- anno[match(rownames(gene_beta), anno$Name), ]
    ord <- order(probe_anno$chr, probe_anno$pos)
    gene_beta <- gene_beta[ord, , drop = FALSE]
    
    # Long format für ggplot
    
    #Long-Format
    df_long <- gene_beta %>%
      as.data.frame() %>%
      rownames_to_column("CpG") %>%
      pivot_longer(cols = -CpG, names_to = "Sample", values_to = "Beta") %>%
      mutate(in_dmr = CpG %in% common_probes)
    
    #Sample_Group hinzufügen
    df_long <- df_long %>%
      left_join(
        dplyr::select(targets, Sample_Name, !!rlang::sym(pd_col)),
        by = c("Sample" = "Sample_Name")
      )
    
    pd_col <- pd_col
    g1 <- compare_values[1]
    g2 <- compare_values[2]
    
    df_long <- df_long %>%
      filter(.data[[pd_col]] %in% compare_values) %>%
      mutate(Compare_Group = factor(.data[[pd_col]], levels = compare_values))
    
    #CpG-Reihenfolge: von unten nach oben = genomisch aufsteigend
    df_long$CpG <- factor(df_long$CpG, levels = rev(rownames(gene_beta)))
    
    label_df <- df_long %>%
      distinct(CpG, in_dmr) %>%
      mutate(
        CpG_label = ifelse(
          in_dmr,
          paste0("<span style='color:red;'>", CpG, "</span>"),
          CpG
        )
      )
    
    sample_order_by_group <- df_long %>%
      group_by(!!rlang::sym(pd_col)) %>%
      group_modify(~ tibble(Sample = get_order_per_group(.x))) %>%
      ungroup()
    
    sample_order <- sample_order_by_group$Sample
    df_long <- df_long %>%
      mutate(Sample = factor(Sample, levels = sample_order))
    
    cpg_levels <- levels(df_long$CpG)
    
    dmr_levels <- cpg_levels[cpg_levels %in% common_probes]
    ymin <- which(cpg_levels == min(dmr_levels))
    ymax <- which(cpg_levels == max(dmr_levels))
    
    p <- ggplot(df_long, aes(x = Sample, y = CpG, fill = Beta)) +
      geom_tile() +
      scale_x_discrete(expand = expansion(add = 0)) +
      scale_fill_gradientn(
        colors = c("yellow", "black", "blue"),
        limits = c(0, 1),
        name = "β-Wert"
      ) +
      facet_grid(cols = vars(.data[[pd_col]]), scales = "free_x", space = "free_x")  +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.title = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      labs(
        x = "Samples",
        y = "CpGs der DMR +- 100 bp",
        title = paste("DMR", i)
      )
    
    if (length(dmr_levels) > 0) {
      ypos <- match(dmr_levels, cpg_levels)
      ymin <- min(ypos) - 0.5
      ymax <- max(ypos) + 0.5
      
      rect_df <- df_long %>%
        distinct(Sample, group = .data[[pd_col]]) %>%   # group ist der Facet-Wert
        count(group, name = "n") %>%
        mutate(
          xmin = 0.5,
          xmax = n + 0.5 - 1,
          ymin = ymin,
          ymax = ymax
        )
      
    }
    
    print(p)
  }
  dev.off()
  
}

makeHeatmap_genelist <- function(pdfname, anno, genelist, myNorm, targets, pd_col, compare_values) {
  pdf(pdfname, width = 10, height = 8)
  i <- 1
  for (gene_of_interest in genelist) {
    
    cat("Starte Heatmap für", gene_of_interest, "\n")
    
    # CpGs des Genes auswählen
    cpgs_for_gene <- anno[grepl(paste0("(^|;)", gene_of_interest, "(;|$)"), anno$UCSC_RefGene_Name, ignore.case=TRUE), ]
    gene_probes_all <- cpgs_for_gene$Name
    common_gene_probes_all <- intersect(gene_probes_all, rownames(myNorm))
    
    if(length(common_gene_probes_all) == 0) {
      cat("Kein CpG für", gene_of_interest, "gefunden – übersprungen.\n")
      next
    }
    
    gene_beta <- myNorm[common_gene_probes_all, , drop = FALSE]
    # CpGs nach genomischer Position sortieren
    
    probe_anno <- anno[match(rownames(gene_beta), anno$Name), ]
    ord <- order(probe_anno$chr, probe_anno$pos)
    gene_beta <- gene_beta[ord, , drop = FALSE]
    # Long format für ggplot
    
    #Long-Format
    df_long <- gene_beta %>%
      as.data.frame() %>%
      rownames_to_column("CpG") %>%
      pivot_longer(cols = -CpG, names_to = "Sample", values_to = "Beta") %>%
      mutate(in_dmr = CpG %in% common_gene_probes_all)
    
    #Sample_Group hinzufügen
    df_long <- df_long %>%
      left_join(
        dplyr::select(targets, Sample_Name, !!rlang::sym(pd_col)),
        by = c("Sample" = "Sample_Name")
      )
    
    pd_col <- pd_col
    g1 <- compare_values[1]
    g2 <- compare_values[2]
    
    df_long <- df_long %>%
      filter(.data[[pd_col]] %in% compare_values) %>%
      mutate(Compare_Group = factor(.data[[pd_col]], levels = compare_values))
    
    #CpG-Reihenfolge: von unten nach oben = genomisch aufsteigend
    df_long$CpG <- factor(df_long$CpG, levels = rev(rownames(gene_beta)))
    
    label_df <- df_long %>%
      distinct(CpG, in_dmr) %>%
      mutate(
        CpG_label = ifelse(
          in_dmr,
          paste0("<span style='color:red;'>", CpG, "</span>"),
          CpG
        )
      )
    
    # Das mit reingenommen:
    sample_order_by_group <- df_long %>%
      group_by(!!rlang::sym(pd_col)) %>%
      group_modify(~ tibble(Sample = get_order_per_group(.x))) %>%
      ungroup()
    
    # Das ersetzt: 
    #Reihenfolge der Samples festlegen
    # sample_order <- df_long %>%
    #   distinct(Sample, .data[[pd_col]]) %>%
    #   arrange(.data[[pd_col]], Sample) %>%
    #   pull(Sample)
    sample_order <- sample_order_by_group$Sample
    
    
    
    # Das ersetzt:
    #df_long$Sample <- factor(df_long$Sample, levels = sample_order)
    df_long <- df_long %>%
      mutate(Sample = factor(Sample, levels = sample_order))
    
    
    cpg_levels <- levels(df_long$CpG)
    
    p <- ggplot(df_long, aes(x = Sample, y = CpG, fill = Beta)) +
      geom_tile() +
      scale_x_discrete(expand = expansion(add = 0)) +
      scale_fill_gradientn(
        colors = c("yellow", "black", "blue"),
        limits = c(0, 1),
        name = "β-Wert"
      ) +
      # Und das ersetzt:
      facet_grid(cols = vars(.data[[pd_col]]), scales = "free_x", space = "free_x")  +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.y = element_text(size = 4),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title = element_text(size = 10),
        panel.grid = element_blank(),
        legend.position = "bottom",
        plot.title = element_text(hjust = 0.5, face = "bold")
      ) +
      labs(
        x = "Samples",
        y = "CpGs der Gene",
        title = paste("Heatmap for Gene:", gene_of_interest)
      )
    
    print(p)
  }
  dev.off()
  
}

get_order_per_group <- function(df) {
  #####
  # help function for ordering the samples in the heatmap-functions based on mean of the beta-values
  #####
  
  mat <- df %>%
    dplyr::select(CpG, Sample, Beta) %>%
    distinct() %>%
    pivot_wider(names_from = Sample, values_from = Beta)
  
  cpg_ids <- mat$CpG
  mat <- as.matrix(mat[,-1, drop = FALSE])
  rownames(mat) <- cpg_ids
  
  if (ncol(mat) <= 1) return(colnames(mat))  # nichts zu clustern
  
  if (anyNA(mat)) {
    row_med <- apply(mat, 1, median, na.rm = TRUE)
    for (r in seq_len(nrow(mat))) {
      na_idx <- is.na(mat[r, ])
      if (any(na_idx)) mat[r, na_idx] <- row_med[r]
    }
  }
  
  hc <- hclust(dist(t(mat)), method = "complete")
  hc$labels[hc$order]
}

annotateDMRs <- function(dmrs, anno, myNorm = NULL) {
  
  #####
  # Annotation von bumphunter-DMRs ausschließlich über CpG-Annotation
  
  # dmrs:
  # Dataframe / bumphunter-Ergebnis mit:
  # - seqnames
  # - start
  # - end
  
  # anno:
  # CpG-Annotation mit mindestens:
  # - Name
  # - chr
  # - pos
  # - UCSC_RefGene_Name
  
  # myNorm:
  # Optional. Matrix mit CpGs als rownames.
  # Falls angegeben, werden nur CpGs berücksichtigt,
  # die auch tatsächlich in myNorm vorhanden sind.
  #####
  
  # bumphunter-Ergebnis in Dataframe umwandeln
  dmrs <- as.data.frame(dmrs)
  
  ########################################
  # Chromosomen vereinheitlichen
  ########################################
  
  dmrs$seqnames <- paste0(
    "chr",
    gsub("^chr", "", as.character(dmrs$seqnames))
  )
  
  anno$chr <- paste0(
    "chr",
    gsub("^chr", "", as.character(anno$chr))
  )
  
  ########################################
  # Neue Annotationsspalten vorbereiten
  ########################################
  
  dmrs$n_CpGs <- 0
  dmrs$CpGs <- NA
  dmrs$Genes <- NA
  
  ########################################
  # DMRs annotieren
  ########################################
  
  for (i in seq_len(nrow(dmrs))) {
    # Aktuelles DMR
    dmr_chr <- dmrs$seqnames[i]
    dmr_start <- dmrs$start[i]
    dmr_end <- dmrs$end[i]
    
    ########################################
    # CpGs innerhalb des DMRs finden
    ########################################
    
    dmr_cpgs <- anno[
      anno$chr == dmr_chr &
        anno$pos >= dmr_start &
        anno$pos <= dmr_end,
      ,
      drop = FALSE
    ]
    
    ########################################
    # Optional:
    # Nur CpGs verwenden, die in myNorm sind
    ########################################
    
    if (!is.null(myNorm)) {
      dmr_cpgs <- dmr_cpgs[
        dmr_cpgs$Name %in% rownames(myNorm),
        ,
        drop = FALSE
      ]
    }
    
    ########################################
    # Keine CpGs gefunden
    ########################################
    
    if (nrow(dmr_cpgs) == 0) {
      cat(
        "Kein CpG für DMR", i,
        "gefunden.\n"
      )
      next 
    }
    ########################################
    # Anzahl CpGs
    ########################################
    
    dmrs$n_CpGs[i] <- nrow(dmr_cpgs)
    
    ########################################
    # CpG-IDs speichern
    ########################################
    
    dmrs$CpGs[i] <- paste(
      unique(dmr_cpgs$Name),
      collapse = ";"
    )
    ########################################
    # Gene aus UCSC_RefGene_Name extrahieren
    ########################################
    
    genes <- unique(
      unlist(
        strsplit(
          dmr_cpgs$UCSC_RefGene_Name,
          ";"
        )
      )
    )
    # Leere / fehlende Gene entfernen
    genes <- genes[
      !is.na(genes) &
        genes != ""
    ]
    ########################################
    # Gene speichern
    ########################################
    if (length(genes) > 0) {
      dmrs$Genes[i] <- paste(
        genes,
        collapse = ";"
      )  
    }
    cat(
      "DMR", i,
      "annotiert:",
      dmrs$n_CpGs[i],
      "CpGs gefunden.\n"
    ) 
  }
  return(dmrs)
  
}
plot_mds_by <- function(dat, samples, group_col, numPos) {
  grp <- samples[[group_col]]
  n_groups <- length(unique(na.omit(grp)))
  if (n_groups <=8){
    colors <- brewer.pal(8, "Dark2")
  }
  else {
    colors <- qualitative_hcl(n_groups, "Dark2")
  }
  mdsPlot_custom(dat, sampGroups = samples[[group_col]], sampNames = NULL, pch = 16, numPositions = numPos, main = paste("Beta MDS 1000 most variable positions grouped by", group_col), pal = colors, legendPos = "bottomright",  legendNCol = 1
  )
}

champ_QC_custom <- function(beta = myLoad$beta,
                            pheno=myLoad$pd$Sample_Group,
                            mdsPlot=TRUE,
                            densityPlot=TRUE,
                            PDFplot=TRUE,
                            Feature.sel="None",
                            resultsDir="./CHAMP_QCimages/",
                            filename_prefix = "",
                            num_genes = 1000)
{
  message("[===========================]")
  message("[<<<<< ChAMP.QC START >>>>>>]")
  message("-----------------------------")
  ### Prepare Checking ###
  if (!file.exists(resultsDir)) dir.create(resultsDir)
  message("champ.QC Results will be saved in ",resultsDir)
  message("[QC plots will be proceed with ",dim(beta)[1], " probes and ",dim(beta)[2], " samples.]\n")
  if(min(beta,na.rm=TRUE)==0)
  {
    beta[beta==0] <- 0.000001
    message("[",length(which(beta==0))," Zeros dectect in your dataset, will be replaced with 0.000001]\n")
  }
  if(ncol(beta)!=length(pheno)) stop("Dimension of DataSet Samples, pheno and name must be the same. Please check your input.")
  message("<< Prepare Data Over. >>")
  
  if(mdsPlot)
  {
    if(PDFplot){
      pdf(paste(resultsDir,paste(filename_prefix, "mdsPlot.pdf", sep = "_"),sep="/"),width=6,height=4)
      mdsPlot(beta,numPositions=num_genes,sampGroups=pheno,colnames(beta))
      dev.off()
    }
    message("<< plot mdsPlot Done. >>\n")
  }
  
  if(densityPlot)
  {
    if(PDFplot){
      pdf(paste(resultsDir,paste(filename_prefix,"densityPlot.pdf", sep="_"),sep="/"),width=6,height=4)
      densityPlot(beta,sampGroups=pheno,main=paste("Density plot of raw data (",nrow(beta)," probes)",sep=""),xlab="Beta")
      dev.off()
    }
    message("<< Plot densityPlot Done. >>\n")
  }
  
  message("[<<<<<< ChAMP.QC END >>>>>>>]")
  message("[===========================]")
  message("[You may want to process champ.norm() next.]\n")
}

mdsPlot_custom <- function(dat, numPositions = numpos, sampNames = NULL,
                           sampGroups = NULL, xlim, ylim, pch = 1,
                           pal = brewer.pal(8, "Dark2"), legendPos = "bottomright",
                           legendNCol, main = NULL) {
  # Check inputs
  if (is(dat, "MethylSet") || is(dat, "RGChannelSet")) {
    b <- getBeta(dat)
  } else if (is(dat, "matrix")) {
    b <- dat
  } else {
    stop("dat must be an 'MethylSet', 'RGChannelSet', or 'matrix'.")
  }
  if (is.null(main)) {
    main <- sprintf(
      "Beta MDS\n%d most variable positions",
      numPositions)
  }
  
  o <- order(rowVars(b), decreasing = TRUE)[seq_len(numPositions)]
  d <- dist(t(b[o, ]))
  fit <- cmdscale(d)
  if (missing(xlim)) xlim <- range(fit[, 1]) * 1.2
  if (missing(ylim)) ylim <- range(fit[, 2]) * 1.2
  if (is.null(sampGroups)) sampGroups <- rep(1, numPositions)
  sampGroups <- as.factor(sampGroups)
  col <- pal[sampGroups]
  if (is.null(sampNames)) {
    plot(
      x = fit[, 1],
      y = fit[, 2],
      col = col,
      pch = pch,
      xlim = xlim,
      ylim = ylim,
      xlab = "",
      ylab = "",
      main = main)
  } else {
    plot(
      x = 0,
      y = 0,
      type = "n",
      xlim = xlim,
      ylim = ylim,
      xlab = "",
      ylab = "",
      main = main)
    text(x = fit[, 1], y = fit[, 2], sampNames, col = col)
  }
  numGroups <- length(levels(sampGroups))
  if (missing(legendNCol)) legendNCol <- numGroups
  if (numGroups > 1) {
    legend(
      x = legendPos,
      cex = 0.6,
      legend = levels(sampGroups),
      ncol = legendNCol,
      text.col = pal[seq_len(numGroups)])
  }
}

SVD_custom <- function(beta=myNorm,
                       rgSet=NULL,
                       pd=myLoad$pd,
                       RGEffect=FALSE,
                       PDFplot=TRUE,
                       resultsDir="./CHAMP_SVDimages/",
                       filename = "SVsummary.pdf")
{
  message("[===========================]")
  message("[<<<<< ChAMP.SVD START >>>>>]")
  message("-----------------------------")
  
  ### Defind some functions to be used in champ.SVD.
  GenPlot_cust <- function(thdens.o,estdens.o,evalues.v){
    minx <- min(min(thdens.o$lambda),min(evalues.v));
    maxx <- max(max(thdens.o$lambda),max(evalues.v));
    miny <- min(min(thdens.o$dens),min(estdens.o$y));
    maxy <- max(max(thdens.o$dens),max(estdens.o$y));
  }
  # Estimate how many latent variable hidden in the data set.
  EstDimRMTv2_cust <- function(data.m)
  {    
    ### standardise matrix
    M <- data.m;
    for(c in 1:ncol(M)) M[,c] <- (data.m[,c]-mean(data.m[,c]))/sqrt(var(data.m[,c]));
    sigma2 <- var(as.vector(M));
    Q <- nrow(data.m)/ncol(data.m);
    thdens.o <- thdens_cust(Q,sigma2,ncol(data.m));
    C <- 1/nrow(M) * t(M) %*% M;
    
    eigen.o <- eigen(C,symmetric=TRUE);
    estdens.o <- density(eigen.o$values,from=min(eigen.o$values),to=max(eigen.o$values),cut=0);
    
    GenPlot_cust(thdens.o,estdens.o,eigen.o$values);
    intdim <- length(which(eigen.o$values > thdens.o$max));
    return(list(cor=C,dim=intdim,estdens=estdens.o,thdens=thdens.o));
  }
  thdens_cust <- function(Q,sigma2,ns)
  {
    lambdaMAX <- sigma2*(1+1/Q + 2*sqrt(1/Q));
    lambdaMIN <- sigma2*(1+1/Q - 2*sqrt(1/Q));
    delta <- lambdaMAX - lambdaMIN;#  print(delta);
    roundN <- 3;
    step <- round(delta/ns,roundN);
    while(step==0){
      roundN <- roundN+1;
      step <- round(delta/ns,roundN);
    }
    lambda.v <- seq(lambdaMIN,lambdaMAX,by=step);
    dens.v <- vector();
    ii <- 1;
    for(i in lambda.v){
      dens.v[ii] <- (Q/(2*pi*sigma2))*sqrt( (lambdaMAX-i)*(i-lambdaMIN) )/i;
      ii <- ii+1;
    }
    return(list(min=lambdaMIN,max=lambdaMAX,step=step,lambda=lambda.v,dens=dens.v));
  }
  
  # A Function to draw heatmap, based on significance level of correlation between phenotype and SVD componnents.
  drawheatmap <- function(svdPV.m) {
    
    op <- par(no.readonly = TRUE)
    on.exit(par(op))
    
    myPalette <- c("darkred","red","orange","pink","white")
    breaks.v  <- c(-10000,-10,-5,-2,log10(0.05),0)
    
    # Mehr Platz: unten für gedrehte Labels, rechts für die Legende
    par(
      mar = c(7, 9, 4, 10) + 0.1,
      mgp = c(3, 1.5, 0),
      xpd = NA
    )
    
    image(
      x = 1:nrow(svdPV.m),
      y = 1:ncol(svdPV.m),
      z = log10(svdPV.m),
      col = myPalette,
      breaks = breaks.v,
      xlab = "", ylab = "",
      axes = FALSE,
      main = "Singular Value Decomposition Analysis (SVD)"
    )
    
    axis(1, at = 1:nrow(svdPV.m),
         labels = paste0("PC-", 1:nrow(svdPV.m)), las = 2)
    
    suppressWarnings(axis(2, at = 1:ncol(svdPV.m),
                          labels = colnames(svdPV.m), las = 2))
    
    usr <- par("usr")
    legend(
      x = usr[2] + 0.5,  # rechts neben den Plot
      y = usr[4],        # oben
      xjust = 0, yjust = 1,
      legend = c(expression("p < 1x"~10^{-10}),
                 expression("p < 1x"~10^{-5}),
                 "p < 0.01", "p < 0.05", "p > 0.05"),
      fill = myPalette,
      bty = "n"
    )
  }
  
  # Function to organize a list for ggplot the screeplot.
  splot_cust <- function(x=svd.o,y=rmt.o)  
  {
    scp <- list(u=x$u[1:nrow(x$u),1:y$dim],v=x$v[1:y$dim,1:y$dim],d=x$d[1:y$dim])
    return(invisible(capture.output(svd.scree_cust(scp))))
  }
  # Function to draw screeplot. Contributed by Rasmus.
  svd.scree_cust <- function(svd.obj, maintitle="Scree Plot", axis.title.x="Component Index (Singular Vectors)", axis.title.y="Percent Variance Explained") 
  {
    if(is.list(svd.obj) & all(names(svd.obj) %in% c("u","d","v"))) {
      print("Your input data is treated as a SVD output, with u, d, v corresponding to left singular vector, singular values, and right singular vectors, respectively.")
    } else {
      print("Your input data is treated as a vector of singular values. For example, it should be svd.obj$d from a SVD output.")
      svd.obj = list(d=svd.obj)
    }
    print("Scree Plot")
    if(is.null(names(svd.obj$d))) 
    {
      names(svd.obj$d) = paste0("Component",1:ncol(svd.obj$v))
    }
    pve = svd.obj$d^2/sum(svd.obj$d^2) * 100
    pve = data.frame(names=names(svd.obj$d), pve)
    pve$names = factor(pve$names, levels=unique(names(svd.obj$d)))
    g = ggplot2::ggplot(pve, aes(names, pve)) + geom_point() + theme_bw()
    gout = g + ylim(0,NA) + 
      ggtitle(maintitle) + # for the main title
      xlab(axis.title.x) + # for the x axis label
      ylab(axis.title.y) + # for the y axis label
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )
    
    return(gout)
  }
  
  ################ Main Program starts here ########################
  
  if (!file.exists(resultsDir)) dir.create(resultsDir)
  message("champ.SVD Results will be saved in ",resultsDir," .\n")
  
  if(length(which(is.na(beta)))>0) message(length(which(is.na(beta)))," NA are detected in your beta Data Set, which may cause fail or uncorrect of SVD analysis. You may want to impute NA with champ.impute() function first.")
  
  if(inherits(beta, "data.frame"))
  {
    message("Your beta parameter is data.frame format. ChAMP is now changing it to matrix.")
    beta <- as.matrix(beta)
  }
  
  message("[SVD analysis will be proceed with ",dim(beta)[1], " probes and ",dim(beta)[2], " samples.]\n")
  
  message("\n[ champ.SVD() will only check the dimensions between data and pd, instead if checking if Sample_Names are correctly matched (because some user may have no Sample_Names in their pd file),thus please make sure your pd file is in accord with your data sets (beta) and (rgSet).]\n")
  
  
  ################ Customise Phenotype Data ########################
  
  if(is.null(pd) | class(pd)=="list") stop("pd parameter in Data Frame or Matrix is necessary And must contain at least tow factors. If your pd is a list, please change its Format.")
  if(class(pd)=="matrix") pd <- as.data.frame(pd)
  
  PhenoTypes.lv_tmp <- pd[,!colnames(pd) %in% c("Sample_Name","Project","filenames","Basename") & apply(pd,2,function(x) length(unique(x)))!=1]
  PhenoTypes.lv <- PhenoTypes.lv_tmp
  
  if(!is.null(rownames(pd))) rownames(PhenoTypes.lv) <- rownames(pd)
  
  if(ncol(PhenoTypes.lv)>=2)
  {
    message("<< Following Factors in your pd(sample_sheet.csv) will be analysised: >>")
    sapply(colnames(PhenoTypes.lv_tmp),function(x) message("<",x,">(",class(PhenoTypes.lv[[x]]),"):",paste(unique(PhenoTypes.lv_tmp[,x]),collapse=", ")))
    message("[champ.SVD have automatically select ALL factors contain at least two different values from your pd(sample_sheet.csv), if you don't want to analysis some of them, please remove them manually from your pd variable then retry champ.SVD().]")
  }else
  {
    stop("You don't have even one factor with at least two value to be analysis. Maybe your factors contains only one value, no variation at all...")
  }
  
  if(ncol(pd) > ncol(PhenoTypes.lv))
  {
    message("\n<< Following Factors in your pd(sample_sheet.csv) will not be analysis: >>")
    sapply(setdiff(colnames(pd),colnames(PhenoTypes.lv)),function(x) message("<",x,">"))
    message("[Factors are ignored because they only indicate Name or Project, or they contain ONLY ONE value across all Samples.]")
  }
  
  #### PhenoTypes.lv prepare ready.
  if(RGEffect==TRUE & is.null(rgSet)) message("If you want to check Effect of Control Probes, you MUST provide rgSet parameter. Now champ.SVD can only analysis factors in pd.")
  if(!is.null(rgSet) & RGEffect)
  {
    if(rgSet@annotation[1]=="IlluminaHumanMethylation450k") data(ControlProbes450K) else data(ControlProbesEPIC)
    dataC2.m <- as.data.frame(log2(apply(ControlProbes,1,function(x) if(x[3]=="Grn") getGreen(rgSet)[x[2],] else getRed(rgSet)[x[2],])))
    PhenoTypes.lv <- cbind(PhenoTypes.lv,dataC2.m)
    message("\n<< Following rgSet information have been added to PhenoTypes.lv. >>")
    sapply(colnames(dataC2.m),function(x) message("<",x,">"))
    #### RG Information extract!
  }	
  
  if(nrow(PhenoTypes.lv)==ncol(beta)) message("\n<< PhenoTypes.lv generated successfully. >>") else stop("Dimension of your pd file (and rgSet information) is not equal to your beta matrix.")
  
  ######################## Do SVD #############################
  
  tmp.m <- beta-rowMeans(beta)
  rmt.o <- EstDimRMTv2_cust(tmp.m);
  svd.o <- svd(tmp.m);
  if(rmt.o$dim > 20) topPCA <- 20  else topPCA <- rmt.o$dim
  
  svdPV.m <- matrix(nrow=topPCA,ncol=ncol(PhenoTypes.lv));
  colnames(svdPV.m) <- colnames(PhenoTypes.lv);
  
  for(c in 1:topPCA)
    for(f in 1:ncol(PhenoTypes.lv))
      if(class(PhenoTypes.lv[,f])!="numeric")
        svdPV.m[c,f] <- kruskal.test(svd.o$v[,c] ~ as.factor(PhenoTypes.lv[[f]]))$p.value
  else
    svdPV.m[c,f] <- summary(lm(svd.o$v[,c] ~ PhenoTypes.lv[[f]]))$coeff[2,4];
  
  message("<< Calculate SVD matrix successfully. >>")
  
  ######################## Plot SVD Image #############################
  
  #Screeplot end
  if(PDFplot)
  {
    pdf(paste(resultsDir,filename,sep=""),width=14,height=10);
    drawheatmap(svdPV.m)
    splot_cust(svd.o,rmt.o)
    dev.off();
  }
  
  message("<< Plot SVD matrix successfully. >>")
  
  message("[<<<<<< ChAMP.SVD END >>>>>>]")
  message("[===========================]")
  message("[If the batch effect is not significant, you may want to process champ.DMP() or champ.DMR() or champ.BlockFinder() next, otherwise, you may want to run champ.runCombat() to eliminat batch effect, then rerun champ.SVD() to check corrected result.]\n")
  return(svdPV.m)
}