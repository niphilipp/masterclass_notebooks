# Nötige Pakete:
library(IlluminaHumanMethylationEPICanno.ilm10b4.hg19) # Annotation for EPICv1 CpGs
library(IlluminaHumanMethylation450kanno.ilmn12.hg19) # Annotation for 450K CpGs
library(dplyr)
library(ChIPseeker)
library(TxDb.Hsapiens.UCSC.hg19.knownGene)
library(tidyr)
library(tibble)

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
    .data[["p.valueArea"]] <= .env$pval,
    .data[["fwerArea"]]    <= .env$fwer,
    abs(.data[["value"]])  >= .env$min_abs_value
  )
  bh_filt %>% 
    arrange(bh_filt$fwerArea)
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
annotateGR <- function(dmrs_gr){
  #####
  # Function for annotating the dmrs calculated by bumphunter
  # dmrs_gr = table / Dataframe of bumphunter results
  #####
  
  dmrs_gr_df <-as.data.frame(dmrs_gr)
  GR <- GRanges(
    seqnames = paste0("chr", gsub("^chr", "", dmrs_gr_df$seqnames)),  # sicher "chr1", "chrX" etc.
    ranges   = IRanges(start = dmrs_gr_df$start, end = dmrs_gr_df$end),
    strand   = "*"
  )
  valid_chr <- paste0("chr", c(1:22, "X", "Y"))
  GR <- GR[seqnames(GR) %in% valid_chr]
  
  GR_anno <- annotatePeak(
    GR,
    tssRegion = c(-2000, 500),          # Promoter = ±2 kb vom TSS
    TxDb      = TxDb.Hsapiens.UCSC.hg19.knownGene,
    annoDb    = "org.Hs.eg.db",
    verbose   = FALSE
  )
  
  GR_prep <- as.data.frame(GR_anno)
  
  dmrs_gr_df$seqnames <- ifelse(grepl("^chr", dmrs_gr_df$seqnames), dmrs_gr_df$seqnames, paste0("chr", dmrs_gr_df$seqnames))
  print(dmrs_gr_df$seqnames)
  GR_prep$seqnames <- as.character(GR_prep$seqnames)
  
  
  # Annotation zu dmr hinzufügen
  dmrs_annot <- dmrs_gr_df %>%
    left_join(
      GR_prep,
      by = c("seqnames", "start", "end")
    )
}

makeHeatmap_dmrregions_cluster <- function(pdfname, dmrs, anno, myNorm, targets, pd_col, compare_values) {
  ##############
  # pdfname = Filename for pdf
  # dmrs = Dataframe of annotated dmrs
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
    dmr_chr <- dmr_row$seqnames
    dmr_start <- dmr_row$start
    dmr_end <- dmr_row$end
    region_start <- dmr_start - 100
    region_end <- dmr_end + 100
    
    cat("Starte Heatmap für DMR", i, "\n")
    
    #CpGs im DMR
    dmr_probes <- anno$Name[
      anno$chr == dmr_chr &
        anno$pos >= dmr_start &
        anno$pos <= dmr_end
    ]
    
    common_probes <- intersect(dmr_probes, rownames(myNorm))
    
    ###########################
    region_probes <- anno$Name[
      anno$chr == dmr_chr &
        anno$pos >= region_start &
        anno$pos <= region_end
    ]
    common_gene_probes_all <- intersect(region_probes, rownames(myNorm))
    
    if(length(common_gene_probes_all) == 0) {
      cat("Kein CpG für DMR ", i, "gefunden – übersprungen.\n")
      next
    }
    
    #Beta-Werte extrahieren & nach Position sortieren
    gene_beta <- myNorm[common_gene_probes_all, , drop = FALSE]
    
    #2️⃣ CpGs nach genomischer Position sortieren
    probe_anno <- anno[match(rownames(gene_beta), anno$Name), ]
    ord <- order(probe_anno$chr, probe_anno$pos)
    gene_beta <- gene_beta[ord, , drop = FALSE]
    
    #3️⃣ Long format für ggplot
    
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
      
      p <- p +
        geom_rect(
          data = rect_df,
          aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = NA,
          color = "darkred",
          linewidth = 0.7
        )
    }
    
    print(p)
  }
  dev.off()
  
}

makeHeatmap_genename_cluster <- function(pdfname, dmrs, anno, myNorm, targets, pd_col, compare_values) {
  ##############
  # pdfname = Filename for pdf
  # dmrs = Dataframe of annotated dmrs
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
    dmr_chr <- dmr_row$seqnames
    dmr_start <- dmr_row$start
    dmr_end <- dmr_row$end
    gene_symbol <- ifelse(is.na(dmr_row$SYMBOL) || dmr_row$SYMBOL == "",
                          paste0("DMR_", i),
                          dmr_row$SYMBOL)
    
    cat("Starte Heatmap für DMR", i, gene_symbol, "\n")
    
    #CpGs im DMR
    dmr_probes <- anno$Name[
      anno$chr == dmr_chr &
        anno$pos >= dmr_start &
        anno$pos <= dmr_end
    ]
    
    common_probes <- intersect(dmr_probes, rownames(myNorm))
    
    ###########################
    cpgs_for_gene <- anno[grepl(paste0("(^|;)", gene_symbol, "(;|$)"), anno$UCSC_RefGene_Name, ignore.case=TRUE), ]
    gene_probes_all <- cpgs_for_gene$Name
    common_gene_probes_all <- intersect(gene_probes_all, rownames(myNorm))
    
    if(length(common_gene_probes_all) == 0) {
      cat("Kein CpG für", gene_symbol, "gefunden – übersprungen.\n")
      next
    }
    #Beta-Werte extrahieren & nach Position sortieren
    gene_beta <- myNorm[common_gene_probes_all, , drop = FALSE]
    
    #2️⃣ CpGs nach genomischer Position sortieren
    probe_anno <- anno[match(rownames(gene_beta), anno$Name), ]
    ord <- order(probe_anno$chr, probe_anno$pos)
    gene_beta <- gene_beta[ord, , drop = FALSE]
    
    #3️⃣ Long format für ggplot
    
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
    
    # Das mit reingenommen:
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
        y = "CpGs der Gene",
        title = paste("DMR", i, "–", gene_symbol)
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
      
      p <- p +
        geom_rect(
          data = rect_df,
          aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
          inherit.aes = FALSE,
          fill = NA,
          color = "darkred",
          linewidth = 0.7
        )
    }
    
    print(p)
  }
  dev.off()
  
}