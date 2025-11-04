
setwd("/Users/vladimir/Downloads")


library(purrr)
library(RColorBrewer)
library(patchwork)
library(dplyr)
library(tidyr)
library(stringr)
library(data.table)
library(ggVennDiagram)
library(ggplot2)
library(forcats)
library(cowplot)

# Set custom theme for all plots
theme_set(theme_classic(base_size = 19))
theme_update(axis.text.x = element_text(colour = "black"), 
             axis.text.y = element_text(colour = "black"),
             plot.title = element_text(hjust = 0.5), 
             plot.subtitle = element_text(hjust = 0.5))

#plot(tpm.sorghum$Sobic.003G234200, tpm.sorghum$Sobic.003G234400)
#cor.test(tpm.sorghum$Sobic.003G234200, tpm.sorghum$Sobic.003G234400)
#!/usr/bin/env Rscript


#===============================================================================
# HELPER FUNCTIONS
#===============================================================================

# Safe correlation function that handles errors gracefully
safe_cor <- possibly(
  function(x, y) {
    out <- suppressWarnings(
      cor.test(x, y, use = "pairwise.complete.obs")
    )
    tibble(
      rho = unname(out$estimate),
      p_value = out$p.value,
      n_complete = out$parameter + 2
    )
  },
  otherwise = tibble(rho = NA_real_, p_value = NA_real_, n_complete = 0)
)

# Function to perform correlation analysis for a region
analyze_region_correlations <- function(region_name, chr, pos_center, win_size = 100000) {
  
  cat("Analyzing region", region_name, "...\n")
  
  # Define genomic window
  lo <- pos_center - win_size
  hi <- pos_center + win_size
  
  # Get reference genes in the window
  ref_genes <- genesPos %>% 
    filter(chr == !!chr,
           start <= hi,
           end >= lo) %>%
    pull(gene)
  
  # Get target genes for this region
  gene_list <- hotspots_sorghumRNA_geneInregions %>%
    filter(region == region_name) %>%
    pull(gene) %>%
    as.character()
  
  # Check which genes exist in expression data
  present_refs <- intersect(ref_genes, names(tpm.sorghum))
  present_genes <- intersect(gene_list, names(tpm.sorghum))
  
  if (length(setdiff(ref_genes, present_refs)) > 0) {
    warning("Missing reference genes: ", 
            paste(setdiff(ref_genes, present_refs), collapse = ", "))
  }
  
  if (length(setdiff(gene_list, present_genes)) > 0) {
    warning("Missing target genes: ", 
            paste(setdiff(gene_list, present_genes), collapse = ", "))
  }
  
  # Ensure numeric data
  tpm_numeric <- tpm.sorghum %>%
    mutate(across(all_of(c(present_refs, present_genes)), ~ as.numeric(.x)))
  
  # Run correlations
  cor_results <- map_dfr(present_refs, function(ref) {
    map_dfr(present_genes, function(g) {
      safe_cor(tpm_numeric[[ref]], tpm_numeric[[g]]) %>%
        mutate(ref_gene = ref, gene = g, .before = rho)
    })
  })
  
  # Generate control correlations
  set.seed(20 + match(region_name, c("A", "B", "C")))
  cset <- sample(colnames(tpm.sorghum), 50)
  
  cor_results_control <- map_dfr(cset, function(ref) {
    map_dfr(present_genes, function(g) {
      safe_cor(tpm_numeric[[ref]], tpm_numeric[[g]]) %>%
        mutate(ref_gene = ref, gene = g, .before = rho)
    })
  })
  
  # Adjust p-values
  cor_results <- cor_results %>%
    group_by(ref_gene) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  cor_results_control <- cor_results_control %>%
    group_by(ref_gene) %>%
    mutate(p_adj = p.adjust(p_value, method = "BH")) %>%
    ungroup()
  
  # Calculate significance cutoffs from control
  eps <- .Machine$double.xmin
  avg_log10_ref <- cor_results_control %>% 
    mutate(p_adj_clamped = ifelse(p_adj == 0, eps, p_adj),
           log10_p_adj = -log10(p_adj_clamped)) %>% 
    group_by(ref_gene) %>% 
    summarise(mean_log10_p = mean(log10_p_adj), .groups = "drop")
  
  cutoffs <- quantile(avg_log10_ref$mean_log10_p, probs = c(0.95, 0.99))
  
  return(list(
    correlations = cor_results,
    cutoffs = cutoffs,
    present_refs = present_refs,
    pos_center = pos_center
  ))
}

# Function to create gene annotation summary
create_gene_summary <- function(ref_genes, pos_center, region_name = "A") {
  
  # Get gene positions
  genesPos2 <- genesPos[genesPos$gene %in% ref_genes, ]
  
  # Convert to data.table for efficient joining
  genesDT <- as.data.table(genesPos2)
  markersDT <- as.data.table(genemarkers)
  
  setkey(genesDT, chr, start, end)
  setkey(markersDT, chr, pos)
  
  # Create start/end columns for markers (SNPs are single positions)
  markersDT[, `:=`(start = pos, end = pos)]
  
  # Perform range join
  result <- foverlaps(markersDT, genesDT, 
                      by.x = c("chr", "start", "end"), 
                      by.y = c("chr", "start", "end"), 
                      nomatch = 0)
  
  # Add missing genes if needed
  missing_genes <- setdiff(ref_genes, unique(result$gene))
  if (length(missing_genes) > 0) {
    for (missing_gene in missing_genes) {
      row_to_add <- genesDT[gene == missing_gene, .(chr, gene, start, end)]
      if (nrow(row_to_add) > 0) {
        result <- rbindlist(list(result, row_to_add), fill = TRUE, use.names = TRUE)
      }
    }
  }
  
  # Add expression status
  result$expressed <- ifelse(result$gene %in% ref_genes, "yes", "no")
  
  # Merge with annotations
  result <- as.data.frame(result)
  result2 <- merge(result, sorghum_annotatitonV5_1, by.x = 2, by.y = 1)
  
  # Add cis-eQTL information
  result2$cis.controlled <- ifelse(result2$gene %in% sorghum.cis$gene, "yes", "no")
  result2$dist.TSS.topSNP <- result2$start - pos_center
  
  # Create gene categories based on impact and cis status
  result3 <- result2 %>% 
    select(gene, marker, description, impact, cis.controlled, expressed, 
           start, dist.TSS.topSNP, `Best-hit-arabi-name`, 
           `Best-hit-arabi-defline`, `Best-hit-rice-defline`) %>%
    group_by(gene) %>% 
    mutate(
      category = case_when(
        region_name == "B" & any(impact == "MODERATE" & cis.controlled == "yes") ~ "Moderate/Cis",
        region_name == "B" & any(impact == "MODERATE" & cis.controlled == "no") ~ "Moderate",
        any(impact == "HIGH" & cis.controlled == "yes") ~ "High/Cis",
        any(impact == "HIGH" & cis.controlled == "no") ~ "High",
        any(cis.controlled == "yes") ~ "Cis",
        TRUE ~ NA_character_
      )
    ) %>% 
    ungroup()
  
  # Summarize at gene level
  gene_summary <- result3 %>% 
    group_by(gene) %>% 
    summarise(
      category = first(category),
      cis.controlled = first(cis.controlled),
      expressed = first(expressed),
      dist.TSS.topSNP = first(dist.TSS.topSNP),
      Best.hit.arabi.name = first(`Best-hit-arabi-name`),
      Best.hit.arabi.defline = first(`Best-hit-arabi-defline`),
      Best.hit.rice.defline = first(`Best-hit-rice-defline`),
      .groups = "drop"
    )
  
  return(gene_summary)
}

# Function to create volcano plot
create_volcano_plot <- function(cor_results, gene_summary, region_name, cutoff, 
                                color_values = NULL) {
  
  # Default colors by region
  if (is.null(color_values)) {
    color_values <- switch(region_name,
                           "A" = c("High/Cis" = "#E41A1C", "High" = "#377EB8", "Cis" = "#4DAF4A"),
                           "B" = c("Moderate/Cis" = "#FF7F00", "Moderate" = "#984EA3", "Cis" = "#4DAF4A"),
                           "C" = c("High/Cis" = "#E41A1C", "High" = "#377EB8", "Cis" = "#4DAF4A")
    )
  }
  
  # Merge correlation data with gene annotations
  plot_df <- cor_results %>% 
    left_join(
      gene_summary %>% select(gene, category, expressed),
      by = c(ref_gene = "gene")
    ) %>% 
    filter(expressed == "yes")
  
  # Simplify gene names for facet labels - keep only last 4 characters
  plot_df$ref_gene <- str_sub(plot_df$ref_gene, -4)
  
  # Calculate mean correlations for each reference gene
  mean_correlations <- plot_df %>%
    mutate(logp = -log10(pmax(p_adj, .Machine$double.xmin))) %>%
    group_by(ref_gene) %>%
    summarise(
      mean_rho = mean(rho, na.rm = TRUE),
      mean_logp = mean(logp, na.rm = TRUE),
      n_points = n(),
      .groups = "drop"
    )
  
  # Debug output
  cat("Mean correlations for region", region_name, ":\n")
  print(mean_correlations)
  
  # Create volcano plot
  p <- plot_df %>% 
    ggplot(aes(x = rho, y = -log10(p_adj), colour = category)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_vline(xintercept = 0, linetype = "dashed") +
    geom_hline(yintercept = cutoff, lty = "dotted", linewidth = 1, color = "black") +
    # Add yellow diamonds for mean correlations
    geom_point(
      data = mean_correlations,
      aes(x = mean_rho, y = mean_logp),
      inherit.aes = FALSE,
      shape = 23,           # filled diamond
      fill = "yellow",
      color = "black",
      size = 3.2,
      stroke = 0.8
    ) +
    facet_wrap(~ ref_gene, nrow = 1, strip.position = "top") +
    scale_x_continuous(
      breaks = c(-1, 0, 1),
      limits = c(-1.05, 1.05),
      expand = expansion(mult = 0.02)
    ) +
    scale_colour_manual(
      values = color_values,
      na.value = "grey70",
      name = "Category"
    ) +
    labs(
      x = "Spearman ρ",
      y = expression(-log[10]("adjusted p"))
    ) +
    theme(
      legend.position = "top",
      strip.placement = "outside",
      strip.text.x = element_text(size = 18, margin = margin(b = 1)),
      strip.background = element_blank()
    )
  
  return(p)
}

#===============================================================================
# DATA LOADING
#===============================================================================

# Load expression data
tpm.sorghum <- fread("counts.NE2021.sorghum.811.filtered.txt", data.table = FALSE)

# Load gene positions
genesPos <- fread("~/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Previous_downloads/Sbicolor_730_v5.1_primary_exons.csv", data.table = FALSE) 
genesPos <- genesPos[genesPos$type == "mRNA", c(1, 3:5)]

# Load hotspot gene assignments
hotspots_sorghumRNA_geneInregions <- fread("hotspots_sorghumRNA_geneInregions.csv", data.table = FALSE)

# Load variant annotations
genemarkers <- fread("SNPeFF.sorghumRNA.csv", data.table = FALSE)
sorghum_annotatitonV5_1 <- fread("sorghum_annotatitonV5_1.csv", data.table = FALSE)

# Load cis-eQTL data
sorghum <- fread("totalcombined_snps_summary_NS_sorghum_2_V2.csv", data.table = FALSE)
sorghum0 <- sorghum[sorghum$num_SNPs >= 3, ]
sorghum.cis <- sorghum0[sorghum0$cis_trans == "cis", ]

#===============================================================================
# REGIONAL ANALYSES
#===============================================================================

# Define hotspot regions
regions <- list(B = list(chr = 10, pos = 53899820))

# Analyze each region
region_results <- map(names(regions), function(region_name) {
  analyze_region_correlations(
    region_name = region_name,
    chr = regions[[region_name]]$chr,
    pos_center = regions[[region_name]]$pos
  )
})
names(region_results) <- names(regions)

# Create gene summaries for each region
gene_summaries <- map2(names(regions), region_results, function(region_name, results) {
  create_gene_summary(
    ref_genes = results$present_refs,
    pos_center = results$pos_center,
    region_name = region_name
  )
})
names(gene_summaries) <- names(regions)

#===============================================================================
# VOLCANO PLOTS
#===============================================================================

# Create volcano plots for each region
volcano_B <- create_volcano_plot(
  region_results$B$correlations,
  gene_summaries$B,
  "B",
  region_results$B$cutoffs[2]
)

#===============================================================================
# DETAILED GENE STRUCTURE PLOTS
#===============================================================================

# Load additional gene structure data
path <- "/Users/vladimir/Library/CloudStorage/OneDrive-UniversityofNebraska-Lincoln/Previous_downloads/"
Sbicolor_730_v5.1_primary_exons <- fread(paste0(path, "Sbicolor_730_v5.1_primary_exons.csv"), data.table = FALSE)

# Standardize UTR names
Sbicolor_730_v5.1_primary_exons$type[Sbicolor_730_v5.1_primary_exons$type == "five_prime_UTR"] <- "5'UTR"
Sbicolor_730_v5.1_primary_exons$type[Sbicolor_730_v5.1_primary_exons$type == "three_prime_UTR"] <- "3'UTR"

# Load trans-eQTL data
sorghum.trans <- sorghum0[sorghum0$cis_trans == "trans", ]
genes_SNPeFF_RNA0.6_topMarker_info <- fread("genes_SNPeFF_RNA0.6_topMarker_info.csv", data.table = FALSE)

# HY5 zoom-in plot (Region B)
create_hy5_zoomin <- function() {
  
  reduced_df4plot <- genes_SNPeFF_RNA0.6_topMarker_info %>%
    filter(region == "B", marker == "Chr10_53899820") %>%
    mutate(chr = chr.x) %>%
    select(marker, chr, pos, genesControled.trans, gene)
  
  reduced_df4plot2 <- sorghum.trans %>%
    filter(CHROM == reduced_df4plot$chr,
           POS >= reduced_df4plot$pos - 10000,
           POS <= reduced_df4plot$pos + 10000)
  
  Sbicolor_730_v5.1_hy5 <- Sbicolor_730_v5.1_primary_exons %>%
    filter(chr == reduced_df4plot$chr,
           start >= reduced_df4plot$pos - 10000,
           end <= reduced_df4plot$pos + 10000,
           type %in% c("5'UTR", "CDS", "3'UTR"))
  
  cisHY5 <- 53899831
  focus_pos_bp <- 53899820
  half_window_kb <- 3
  half_window_bp <- half_window_kb * 1000
  
  ggplot(reduced_df4plot2, aes(x = POS/1e6, y = 200)) +
    geom_rect(
      data = Sbicolor_730_v5.1_hy5,
      inherit.aes = FALSE,
      aes(xmin = start/1e6, xmax = end/1e6,
          ymin = 0, ymax = nrow(reduced_df4plot2) * 0.2,
          fill = type),
      color = "black"
    ) +
    geom_segment(
      data = Sbicolor_730_v5.1_hy5 %>%
        group_by(gene) %>%
        arrange(start) %>%
        mutate(next_start = lead(start), next_end = lead(end)) %>%
        filter(!is.na(next_start)),
      inherit.aes = FALSE,
      aes(x = end/1e6, xend = next_start/1e6,
          y = (nrow(reduced_df4plot2) * 0.2)/2,
          yend = (nrow(reduced_df4plot2) * 0.2)/2),
      color = "black", linewidth = 3
    ) +
    geom_jitter(shape = 6, height = nrow(reduced_df4plot2)/3, width = 0, size = 5, alpha = 0.5) +
    annotate("point", x = cisHY5/1e6, y = 60, shape = 17, size = 5, stroke = 1.1, color = "yellow") +
    scale_x_continuous(
      limits = c((focus_pos_bp - half_window_bp)/1e6, (focus_pos_bp + half_window_bp)/1e6),
      expand = c(0, 0),
      breaks = round(seq(from = (focus_pos_bp - half_window_bp)/1e6,
                         to = (focus_pos_bp + half_window_bp)/1e6,
                         by = 0.005), 3)
    ) +
    scale_fill_manual(values = c("5'UTR" = "skyblue", "3'UTR" = "grey", "CDS" = "blue")) +
    coord_cartesian(xlim = c((focus_pos_bp - half_window_bp)/1e6, (focus_pos_bp + half_window_bp)/1e6)) +
    labs(x = "Genomic Position (Mb)", y = "", fill = "Feature") +
    theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none")
}

# Create zoom-in plots
hy5_zoomin <- create_hy5_zoomin()

###########################################
###########################################
# Libraries
###########################################
suppressPackageStartupMessages({
  library(ape)
  library(pegas)
  library(igraph)
  library(RColorBrewer)
  library(dplyr)
  library(ggplot2)
  library(scatterpie)  # for pie nodes
  library(svglite)     # optional for clean SVG export
})

###########################################
# Inputs and layout knobs
###########################################
args <- commandArgs(trailingOnly = TRUE)
fasta_file    <- "HaplotypeNet/Sobic010G183600_haplotypes_noamb.fasta"
metadata_file <- if (length(args) >= 1) args[1] else "HaplotypeNet//Sobic010G183600_allele_chr10_53899820_group.tsv"

length_power <- 1.0   # amplify long edges (>1) or keep proportional (=1)
kk_const     <- 1.6   # spring constant for KK layout
spread       <- 3.0   # scale coordinates after layout to add whitespace

# ---- Pie sizing knobs ----
R_BASE_FRAC  <- 0.08  # base radius as fraction of median MST edge length
R_SCALE_FRAC <- 0.05  # extra radius scaled by sqrt(freq) range
PIE_SCALE    <- 4   # <--- Global multiplier (like your SHRINK, but for pies). >1 bigger; <1 smaller.

###########################################
# Load FASTA
###########################################
dna <- read.dna(fasta_file, format = "fasta", as.matrix = TRUE)
if (is.null(rownames(dna))) stop("FASTA must contain sample names in headers.")
sample_names <- rownames(dna)

###########################################
# Collapse to haplotypes (ambiguity-tolerant)
###########################################
haps <- haplotype(dna, strict = FALSE)
idx <- attr(haps, "index")  # list: per haplotype samples (names or indices)

hap_names <- rownames(haps)
if (is.null(hap_names)) hap_names <- paste0("Hap", seq_along(idx))

to_samples <- function(x) {
  if (is.numeric(x)) sample_names[x] else unname(as.character(x))
}

freq <- vapply(idx, function(x) length(to_samples(x)), integer(1))

###########################################
# Metadata (sample -> TT/CC/CT) and filtering
###########################################
metadata <- read.table(metadata_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!all(c("sample", "group") %in% names(metadata))) {
  stop("Metadata file must contain columns named 'sample' and 'group'.")
}
group_map <- setNames(metadata$group, metadata$sample)

cats_keep <- c("TT", "CC", "CT")
keep_idx <- vapply(seq_along(idx), function(i) {
  samp <- to_samples(idx[[i]])
  grp  <- group_map[samp]
  sum(!is.na(grp) & grp %in% cats_keep) > 2
}, logical(1))

indices_plot    <- idx[keep_idx]
hap_names_plot  <- hap_names[keep_idx]
freq_plot       <- vapply(indices_plot, function(x) length(to_samples(x)), integer(1))
haps_plot       <- haps[keep_idx, , drop = FALSE]

###########################################
# Distances (Hamming ignoring Ns) & MST
###########################################
dist_mat <- as.matrix(dist.dna(haps_plot, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE))
full_graph <- graph_from_adjacency_matrix(dist_mat, mode = "undirected", weighted = TRUE, diag = FALSE)
mst_graph  <- mst(full_graph, weights = E(full_graph)$weight)

###########################################
# Layout (KK with transformed edge lengths) + radial spacing
###########################################
set.seed(42)
w <- E(mst_graph)$weight
offset_len <- 1.0
L <- offset_len + log1p(w)                 # compress dynamic range slightly
w_layout <- (L + 1e-6) ^ length_power

layout_raw <- layout_with_kk(mst_graph, weights = w_layout, kkconst = kk_const)

sx <- max(layout_raw[,1]) - min(layout_raw[,1]); if (sx == 0) sx <- 1
sy <- max(layout_raw[,2]) - min(layout_raw[,2]); if (sy == 0) sy <- 1
layout_coords <- cbind(
  (layout_raw[,1] - mean(layout_raw[,1]))/sx,
  (layout_raw[,2] - mean(layout_raw[,2]))/sy
) * spread

# anchor largest haplotype to center & radial spacing by MST depth
center_idx <- which.max(freq_plot)
center_xy  <- layout_coords[center_idx, , drop = FALSE]
layout_coords[,1] <- layout_coords[,1] - center_xy[1,1]
layout_coords[,2] <- layout_coords[,2] - center_xy[1,2]

depths <- distances(mst_graph, v = center_idx)
depths <- as.numeric(depths[1, ])
depths[!is.finite(depths)] <- max(depths[is.finite(depths)], na.rm = TRUE)
maxd <- if (all(!is.finite(depths))) 1 else max(depths, na.rm = TRUE)
if (maxd < 1) maxd <- 1
radial_scale <- 1 + 0.25 * (depths / maxd)
for (i in seq_len(nrow(layout_coords))) {
  r <- sqrt(layout_coords[i,1]^2 + layout_coords[i,2]^2)
  if (r > 0) {
    ux <- layout_coords[i,1] / r; uy <- layout_coords[i,2] / r
    layout_coords[i,1] <- ux * r * radial_scale[i]
    layout_coords[i,2] <- uy * r * radial_scale[i]
  }
}

xr <- range(layout_coords[,1]); yr <- range(layout_coords[,2])
padx <- diff(xr)*0.05; pady <- diff(yr)*0.05
xlim <- c(xr[1]-padx, xr[2]+padx); ylim <- c(yr[1]-pady, yr[2]+pady)

###########################################
# Aesthetics / edges + labels
###########################################
edge_widths  <- 1.2 + 2.5 / E(mst_graph)$weight
edge_labels  <- round(E(mst_graph)$weight)
label_lines  <- paste0(hap_names_plot, "\n(", freq_plot, ")")

###########################################
# Build node/edge data frames
###########################################
nodes_df <- tibble(
  id      = seq_len(vcount(mst_graph)),
  x       = layout_coords[,1],
  y       = layout_coords[,2],
  label   = label_lines,
  freq    = freq_plot
)

E_mat <- igraph::as_edgelist(mst_graph, names = FALSE)
edges_df <- tibble(
  from_id = E_mat[,1],
  to_id   = E_mat[,2],
  weight  = E(mst_graph)$weight
) %>%
  left_join(nodes_df %>% select(id, x, y), by = c("from_id" = "id")) %>%
  rename(x = x, y = y) %>%
  left_join(nodes_df %>% select(id, x, y), by = c("to_id" = "id"), suffix = c("", ".to")) %>%
  rename(xend = x.to, yend = y.to) %>%
  mutate(mx = (x + xend)/2,
         my = (y + yend)/2,
         lw = edge_widths,
         edge_lab = as.character(round(weight)))

# median edge length in plot units
edge_len <- sqrt((edges_df$x - edges_df$xend)^2 + (edges_df$y - edges_df$yend)^2)
med_edge <- stats::median(edge_len, na.rm = TRUE)
if (!is.finite(med_edge) || med_edge == 0) med_edge <- 1

# base radius scaled by median edge length; then apply a global multiplier (PIE_SCALE)
sqrtf <- sqrt(nodes_df$freq)
sqrtf_scaled <- if (max(sqrtf) > 0) sqrtf / max(sqrtf) else sqrtf
nodes_df <- nodes_df %>%
  mutate(radius = med_edge * (R_BASE_FRAC + R_SCALE_FRAC * sqrtf_scaled) * PIE_SCALE)

###########################################
# Per-node TT/CC/CT counts for pies
###########################################
get_counts <- function(sample_indices) {
  s <- to_samples(sample_indices)
  g <- group_map[s]
  g <- g[!is.na(g)]
  c(TT = sum(g == "TT"), CC = sum(g == "CC"), CT = sum(g == "CT"))
}
counts_mat <- t(vapply(indices_plot, get_counts, numeric(3)))
colnames(counts_mat) <- c("TT", "CC", "CT")

pie_df <- nodes_df %>%
  bind_cols(as.data.frame(counts_mat)) %>%
  mutate(total = TT + CC + CT)
# If needed, drop empty pies:
# pie_df <- dplyr::filter(pie_df, total > 0)

###########################################
# Colors for slices
###########################################
slice_palette <- c(TT = "#1b9e77", CC = "#d95f02", CT = "#7570b3")

###########################################
# GGPlot with pie nodes
###########################################
if (!exists("xlim") || !exists("ylim")) {
  pad <- 0.05
  xr <- range(nodes_df$x); yr <- range(nodes_df$y)
  xd <- diff(xr); yd <- diff(yr)
  xlim <- xr + c(-pad*xd, pad*xd)
  ylim <- yr + c(-pad*yd, pad*yd)
}

p_net_pies <- ggplot() +
  # edges
  geom_segment(
    data = edges_df,
    aes(x = x, y = y, xend = xend, yend = yend, linewidth = lw),
    color = "black", lineend = "round", alpha = 0.65
  ) +
  # edge labels
  geom_text(data = edges_df, aes(x = mx, y = my, label = edge_lab),
            size = 3.0, color = "black") +
  # pie nodes
  scatterpie::geom_scatterpie(
    data = pie_df,
    aes(x = x, y = y, r = radius),
    cols = c("TT", "CC", "CT"),
    color = "grey30", alpha = 0.98
  ) +
  scale_fill_manual(values = slice_palette, breaks = c("TT","CC","CT"),
                    name = "Allele at Chr10:53899820", drop = FALSE) +
  # node labels on top
  geom_text(
    data = nodes_df,
    aes(x = x, y = y, label = label),
    color = "black", size = 2.7,
    vjust = 0.5, hjust = 0.5
  ) +
  scale_linewidth_identity(guide = "none") +
  coord_equal(xlim = xlim, ylim = ylim, expand = FALSE) +
  # (Optional) add a discrete legend swatch
  guides(fill = guide_legend(title = "Allele at Chr10:53899820")) +
  theme_void(base_family = "Helvetica") +
  theme(
    legend.position = c(0.68, 0.95),
    legend.justification = c(0, 1)
  )

print(p_net_pies)


###############################################
###############################################
###############################################
###############################################

fasta_file <- "HaplotypeNet/Sobic010G183600_haplotypes_noamb.fasta"
allele_tsv <- "HaplotypeNet/Sobic010G183600_allele_chr10_53899820_group.tsv"
outprefix  <- "Sobic010G183600_pairwise"

# Read sequences
dna <- read.dna(fasta_file, format = "fasta", as.matrix = TRUE)
if (is.null(rownames(dna))) stop("FASTA must contain sample names as headers.")
samples <- rownames(dna)

# Sample -> group
a <- read.table(allele_tsv, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
if (!all(c("sample", "group") %in% names(a))) stop("Allele TSV must have columns: sample, group")
a$sample <- trimws(a$sample); a$group <- trimws(a$group)
gmap <- setNames(a$group, a$sample)
grp <- unname(gmap[samples])
grp[is.na(grp) | grp == ""] <- "missing"

# Pairwise distances ignoring Ns
dm <- as.matrix(dist.dna(dna, model = "N", pairwise.deletion = TRUE, as.matrix = TRUE))
if (!all(rownames(dm) == samples)) dm <- dm[samples, samples]

# Long format for i<j
n <- length(samples)
idx_i <- integer(); idx_j <- integer(); dval <- numeric()
for (i in seq_len(n - 1)) {
  jseq <- (i + 1):n
  idx_i <- c(idx_i, rep(i, length(jseq)))
  idx_j <- c(idx_j, jseq)
  dval  <- c(dval, dm[i, jseq])
}
df <- data.frame(
  sample1 = samples[idx_i],
  sample2 = samples[idx_j],
  group1 = grp[idx_i],
  group2 = grp[idx_j],
  dist = dval,
  stringsAsFactors = FALSE
)

# Category labels
cat <- character(nrow(df))
same <- df$group1 == df$group2
cat[same & df$group1 == "CC"] <- "CC"
cat[same & df$group1 == "TT"] <- "TT"
cat[!same & ((df$group1 == "CC" & df$group2 == "TT") | (df$group1 == "TT" & df$group2 == "CC"))] <- "CC/TT"
cat[cat == ""] <- "Other"
df$category <- factor(cat, levels = c("CC", "TT", "CC/TT", "Other"))


df2 <- df %>%
  filter(category != "Other")

# Write long TSV
outfile_tsv <- paste0(outprefix, "_distance_long.tsv")
write.table(df2, file = outfile_tsv, sep = "\t", quote = FALSE, row.names = FALSE)


head(df2)

# Desired order and palette (note the double backslash in strings)
lvl_order <- c("CC", "CC\\TT", "TT")
slice_palette <- c(CC = "#d95f02", `CC\\TT` = "#7570b3", TT = "#1b9e77")

# Robust category mapping using group1/group2, then force order
dfp <- df2 %>%
  mutate(
    category_norm = case_when(
      group1 == "CC" & group2 == "CC" ~ "CC",
      group1 == "TT" & group2 == "TT" ~ "TT",
      (group1 == "CC" & group2 == "TT") | (group1 == "TT" & group2 == "CC") ~ "CC\\TT",
      TRUE ~ "CC\\TT"  # fallback if anything else sneaks in
    ),
    category_norm = factor(category_norm, levels = lvl_order)
  )

# Counts and y-positions per category (use each group's own max for safer label placement)
counts <- dfp %>%
  group_by(category_norm) %>%
  summarize(n = n(), y = max(dist, na.rm = TRUE) * 1.05, .groups = "drop")

# Boxplot
box <- ggplot(dfp, aes(x = category_norm, y = dist, fill = category_norm)) +
  geom_boxplot(outlier.shape = 21, alpha = 0.85, width = 0.6) +
  geom_text(data = counts, aes(y = y, label = paste0(n)), vjust = 0, size = 6) +
  scale_fill_manual(values = slice_palette, drop = FALSE) +
  labs(x = "Category", y = "Pairwise Distance") +
  theme(legend.position = "none")


###############################################
###############################################
###############################################
###############################################

# ==============================================================================
# REGION B ANALYSIS - HY5 VALIDATION (EXACT COPY FROM ORIGINAL)
# ==============================================================================

print("=== REGION B ANALYSIS - HY5 VALIDATION ===")

# Load required data files
hotspots_sorghumRNA_geneInregions <- fread("GeneExplanation/hotspots_sorghumRNA_geneInregions.csv", data.table = F)
sorghum_arabidopsis_gene_function <- fread("Sorghum_arabidopsis_gene_function.csv", data.table = F)
HY5.target <- fread("HY5_Target_Genes.csv", data.table = F)
tpm.sorghum <- fread("GeneExplanation/counts.NE2021.sorghum.811.filtered.txt", data.table = F)

# Calculate overlaps
# Clean arabidopsis gene names
HY5.target$Gene_At <- gsub("At","AT", HY5.target$Gene_At)
HY5.target$Gene_At <- gsub("g","G", HY5.target$Gene_At)

nrow(HY5.target)

SorghumHitsToTAIR <- fread("SorghumHitsToTAIR.tsv", data.table = F)

# V1: "Sobic.K022032.1.p" -> "Sobic.K022032"
SorghumHitsToTAIR[[1]] <- sub('^([^.]+\\.[^.]+).*$', '\\1', as.character(SorghumHitsToTAIR[[1]]))

# V2: "ATMG00990.1" -> "ATMG00990"
SorghumHitsToTAIR[[2]] <- sub('\\..*$', '', as.character(SorghumHitsToTAIR[[2]]))
SorghumHitsToTAIR <- SorghumHitsToTAIR[,c(2,1,3:4)]
colnames(SorghumHitsToTAIR)[2] <- "sorghum"
colnames(SorghumHitsToTAIR)[1] <- "arabidopsis"
colnames(SorghumHitsToTAIR)[3] <- "scores"
colnames(SorghumHitsToTAIR)[4] <- "e_values"

SorghumHitsToTAIR2 <- na.omit(SorghumHitsToTAIR)

nrow(SorghumHitsToTAIR)

head(SorghumHitsToTAIR2)

# keep ONE row per pair, preferring highest score, then lowest e-value
SorghumHitsToTAIR3 <- SorghumHitsToTAIR2 %>%
  arrange(arabidopsis, sorghum, desc(scores), e_values) %>%  # tie-break order
  distinct(arabidopsis, sorghum, .keep_all = TRUE)


head(SorghumHitsToTAIR3)

# 1) For each Arabidopsis gene, how many distinct sorghum genes match?
per_arab_counts <- SorghumHitsToTAIR3 %>%
  distinct(arabidopsis, sorghum) %>%            # avoid double-counting same pair
  count(arabidopsis, name = "n_sorghum_hits") %>%
  arrange(desc(n_sorghum_hits))

# 2) Subset: Arabidopsis genes that match multiple sorghum genes
multi_map <- per_arab_counts %>%
  filter(n_sorghum_hits > 1)

# How many Arabidopsis genes have multiple sorghum matches?
n_multi <- nrow(multi_map)
n_multi

# 3) (Optional) See the actual sorghum IDs per Arabidopsis gene
multi_map_with_ids <- SorghumHitsToTAIR3 %>%
  distinct(arabidopsis, sorghum) %>%
  group_by(arabidopsis) %>%
  summarise(
    n_sorghum_hits = n_distinct(sorghum),
    sorghum_ids = paste(sort(unique(sorghum)), collapse = ";")
  ) %>%
  ungroup() %>%
  filter(n_sorghum_hits > 1) %>%
  arrange(desc(n_sorghum_hits), arabidopsis)

# Peek
head(multi_map_with_ids, 10)

# 4) Distribution of how many sorghum genes per Arabidopsis gene
hit_distribution <- per_arab_counts %>%
  count(n_sorghum_hits, name = "n_arabidopsis_genes") %>%
  arrange(n_sorghum_hits)

hit_distribution

# 5) Example: all sorghum matches for a specific Arabidopsis gene (e.g., AT1G01050)
SorghumHitsToTAIR3 %>%
  filter(arabidopsis == "AT1G01050") %>%
  distinct(arabidopsis, sorghum) %>%
  pull(sorghum)


### expressed in TPM
SorghumHitsToTAIR4 <- SorghumHitsToTAIR3[SorghumHitsToTAIR3$sorghum %in% colnames(tpm.sorghum),]

#now that we have all arabidopsis genes that can match more than one sorghum gene we add a "yes" if there are targets
SorghumHitsToTAIR4$HY5.target <- ifelse(SorghumHitsToTAIR4$arabidopsis %in% HY5.target$Gene_At, "yes", NA)

NROW(SorghumHitsToTAIR4[duplicated(SorghumHitsToTAIR4$sorghum),])

#remove sorghum duplicated genes, now we should have all the information 
SorghumHitsToTAIR5 <- SorghumHitsToTAIR4 %>%
  distinct(sorghum, .keep_all = T)

# Filter hotspots for region B
hotspots_sorghumRNA_geneInregions.B <- hotspots_sorghumRNA_geneInregions %>%
  filter(region == "B")
nrow(hotspots_sorghumRNA_geneInregions.B)

colnames(hotspots_sorghumRNA_geneInregions.B)[9] <- "sorghum"

R.B.ara <- hotspots_sorghumRNA_geneInregions.B %>%
  left_join(SorghumHitsToTAIR5, by = "sorghum") %>%
  filter(!is.na(arabidopsis)) %>%
  distinct(sorghum, .keep_all = TRUE)

nrow(R.B.ara)
337/356

overlap_hy5_B <- nrow(R.B.ara[which(R.B.ara$HY5.target == "yes"),])
hotspot_genes_B <- nrow(R.B.ara)
total_hy5_targets <- nrow(SorghumHitsToTAIR5[which(SorghumHitsToTAIR5$HY5.target == "yes"),])
universe <- nrow(SorghumHitsToTAIR5) # Total genes with arabidopsis annotations

18506/(ncol(tpm.sorghum)-1)

18506-3522

print(paste("HY5 targets in region B:", overlap_hy5_B))
print(paste("Total HY5 targets in genome:", total_hy5_targets))



# Build Fisher's exact test contingency table
contingency.hy5.B <- matrix(
  c(overlap_hy5_B,                                    # in Hotspot & in HY5 targets
    hotspot_genes_B - overlap_hy5_B,                  # in Hotspot & NOT in HY5 targets
    total_hy5_targets - overlap_hy5_B,                # NOT in Hotspot & in HY5 targets: total_hy5_targets - overlap_hy5_B, 
    universe - hotspot_genes_B - total_hy5_targets + overlap_hy5_B),  # neither: universe - hotspot_genes_B - total_hy5_targets + overlap_hy5_B
  nrow = 2,
  byrow = TRUE,
  dimnames = list(
    Hotspot = c("yes","no"),
    HY5_target = c("yes","no"))
)

print("Contingency table for HY5 enrichment:")
print(contingency.hy5.B)

118+219+3404+14765

# Perform Fisher's exact test
ft.hy5.B <- fisher.test(contingency.hy5.B, alternative = "greater")
print("Fisher's exact test results:")
print(ft.hy5.B)

# Create Fisher plot
plot_df <- as.data.frame(as.table(contingency.hy5.B)) %>% 
  dplyr::rename(Count = Freq)  

fisherplot.hy5.B <- ggplot(plot_df, aes(x = Hotspot, y = Count, fill = HY5_target)) +
  geom_col(position = "fill", width = 1, colour = "black") +
  geom_text(aes(label = Count),
            position = position_fill(vjust = 0.5),
            colour   = "white", size = 6) +
  # Add p-value annotation inside the plot
  annotate("text", x = 1.5, y = 0.9, 
           label = paste0("p = ", signif(ft.hy5.B$p.value, 3)), 
           size = 6) +
  scale_fill_manual(
    values = c("yes" = "#FF6B6B", "no" = "#D3D3D3"),  
    labels = c("HY5 targets", "Background")
  ) +
  scale_x_discrete(labels = c("B targets", "Genome")) +
  guides(fill = guide_legend(title = NULL)) +
  labs(y = "Proportion of genes", x = "Category") +
  theme(legend.position = "top",  
        legend.box.margin = margin(b = -30))  # Reduced top margin to bring legend closer
fisherplot.hy5.B



###############################################
###############################################
###############################################
###############################################

results_phewas_ravi_sig <- fread("GeneExplanation/results_cohen_sig2.csv", data.table = F)
ravi_234_categories <- fread("ravi_234_categories.csv", data.table = F)

# Filter phenotype data for topB marker
results_phewas_ravi_sig_topB <- results_phewas_ravi_sig %>%
  filter(marker == "topB")

print(paste("Number of phenotypes significantly associated with topB:", nrow(results_phewas_ravi_sig_topB)))

# Merge with phenotype categories
results_phewas_ravi_sig_topB <- merge(results_phewas_ravi_sig_topB, ravi_234_categories, by = "phenotype")


fwrite(results_phewas_ravi_sig_topB, "results_phewas_ravi_sig_topB.csv")

######### calculate Cohen's D (alt - ref) for the top candidate gene 
results_cohen_sig_geneExp <- fread("results_cohen_sig_geneExp.csv", data.table = F)

head(results_cohen_sig_geneExp)
results_cohen_sig_geneExp$phenotype2 <- c("SbCA2","SbCA2","SbCA2","SbCA2","SbHY5a","SbHY5a","SbHY5a","SbHY5a","SbPAA1-like","SbPAA1-like","SbPAA1-like")
results_cohen_sig_geneExp$category <- c("SbCA2","SbCA2","SbCA2","SbCA2","SbHY5a","SbHY5a","SbHY5a","SbHY5a","SbPAA1-like","SbPAA1-like","SbPAA1-like")


results_cohen_sig_geneExp <- results_cohen_sig_geneExp %>%
  filter(marker == "topB", phenotype == "Sobic.010G183600")

# 1) Define the template (columns you want to end up with)
template_cols <- names(results_phewas_ravi_sig_topB)

# 2) Add any missing columns to results_cohen_sig_geneExp, filled with NA
missing_cols <- setdiff(template_cols, names(results_cohen_sig_geneExp))
if (length(missing_cols)) {
  results_cohen_sig_geneExp[missing_cols] <- NA
}

# 3) Keep only the template columns and in the same order
results_cohen_sig_geneExp_aligned <- results_cohen_sig_geneExp %>%
  select(all_of(template_cols))

# 4) Row-bind
combined <- bind_rows(
  results_phewas_ravi_sig_topB %>% select(all_of(template_cols)),
  results_cohen_sig_geneExp_aligned
)

# combined now has the same schema as results_phewas_ravi_sig_topB

combined
fwrite(combined, "PheWAS_B.csv")


# Create category-level summary plot for topB only
if (nrow(combined) > 0) {
  # Set up for category analysis
  metric_col <- "cohen_d_alt_minus_ref"
  metric_lab <- "Cohen's D (alt - ref)"
  
  # Keep only rows with numeric Cohen's D
  sig <- combined %>% 
    filter(is.finite(.data[[metric_col]])) %>%
    mutate(Block = category)
  
  # Remove rows with missing categories
  sig <- sig[!is.na(sig$Block), ]
  
  if (nrow(sig) > 0) {
    # Compute per-block summary
    stats <- sig %>%
      group_by(Block) %>%
      summarise(
        mean   = mean(.data[[metric_col]], na.rm = TRUE),
        sd     = sd(.data[[metric_col]],   na.rm = TRUE),
        count  = dplyr::n(),
        stderr = sd / sqrt(pmax(count, 1)),
        .groups = "drop"
      )
    
    # Order blocks by mean Cohen's D
    stats <- stats %>% mutate(Block = fct_reorder(Block, mean))
    
    # Data for individual points
    dots <- sig %>%
      transmute(
        Block,
        value = .data[[metric_col]]
      ) %>%
      mutate(Block = factor(Block, levels = levels(stats$Block)))
    
    # Category-level barplot with error bars
    hy5.category.plot <- ggplot(stats, aes(x = Block, y = mean)) +
      geom_col(fill = "gold2", width = 0.7, color = "black", linewidth = 0.2) +
      geom_errorbar(aes(ymin = mean - stderr, ymax = mean + stderr), width = 0.25) +
      geom_hline(yintercept = 0, color = "black") +
      geom_hline(yintercept = c(-0.8, -0.5, -0.2, 0.2, 0.5, 0.8),
                 linetype = "dashed", linewidth = 0.25, alpha = 0.4) +
      geom_point(
        data = dots,
        aes(x = Block, y = value),
        inherit.aes = FALSE,
        shape = 21,              # allows separate fill and outline
        size = 5,
        alpha = 0.8,
        fill = "gold4",           # fill color for circles
        color = "black",         # outline color (optional, can remove if not needed)
        position = position_jitter(width = 0.15, height = 0)
      ) +
      coord_flip() +
      labs(x = NULL, y = metric_lab)
  }
}

hy5.category.plot


#===============================================================================
# FINAL FIGURE ASSEMBLY
#===============================================================================
#volcano_B
#hy5_zoomin
#p_net_pies
#box
#fisherplot.hy5.B
#hy5.category.plot

plots <-
  (volcano_B / hy5_zoomin / box) +
  plot_layout(heights = c(1, 1, 1))

bottom_row <- plot_grid(plots[[2]], p_net_pies, rel_widths = c(1,1),
                        nrow = 1, labels = c('B', 'C'),hjust = 0, vjust = 1, label_size = 18)
bottom_row2 <- plot_grid(plots[[3]],fisherplot.hy5.B, hy5.category.plot, rel_widths = c(.8, .8,1),
                         nrow = 1, labels = c('D','E', 'F'),hjust = 0, vjust = 1, label_size = 18)

fancyplots2 <- plot_grid(plots[[1]], bottom_row, bottom_row2, ncol = 1,  labels = c('A', ''),hjust = 0, vjust = 1, label_size = 18)
fancyplots2


ggsave(plot=fancyplots2, "SbHY5a_plot.svg", width = 12, height = 11)

