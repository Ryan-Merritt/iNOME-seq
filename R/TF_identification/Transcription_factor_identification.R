###############
# Preparation #
###############

# Setting working directory
setwd("D:/iNOME/")

# Loading required libraries
library("cowplot")
library("circlize")
library("lemon")
library("ggforce")
library("ggrepel")
library("patchwork")
library("rtracklayer")
library("DMRcaller")
library("patchwork")
library("ChIPpeakAnno")
library("EnrichedHeatmap")
library("ggpubr")
library("GenomicFeatures")
library("tidyverse")
library("plyranges")
library("slider")
library("VennDiagram")
library("regioneR")
library("ggplotify")
library("drfun")

# Choosing tissue
tissue <- "root"

# Choosing pipeline
pipeline <- "merged" # "merged" "split"
pipeline_id <- "M" # "M" "S"

######
# GC #
######

# Loading reference genome
genome <- import_genome("reference/TAIR10_chr_all.fasta")

# Loading CX reports
if (tissue == "leaf") {
    load("data/cx_report/leaf/leafMethylationData.RData")
    control <- leafMetData$CE_24h_L.R1.L6.
    control <- dropSeqlevels(control, c("Pt", "Mt"), "coarse")
    seqlevels(control) <- as.character(1:5)

    mcvip <- leafMetData$"ME_24h_L.R1.L6."
    mcvip <- dropSeqlevels(mcvip, c("Pt", "Mt"), "coarse")
    seqlevels(mcvip) <- as.character(1:5)
} else {
    load("data/cx_report/root/rootMethylationData.RData")
    control <- rootMetData$CE_12h_R.R1.RUN0076.L3.
    control <- dropSeqlevels(control, c("Pt", "Mt"), "coarse")
    seqlevels(control) <- as.character(1:5)

    mcvip <- rootMetData$"ME_12h_R.R1.RUN0076.L3."
    mcvip <- dropSeqlevels(mcvip, c("Pt", "Mt"), "coarse")
    seqlevels(mcvip) <- as.character(1:5)
}

# Extracting GC
control_gc <- extract_sequence(control, genome)
mcvip_gc <- extract_sequence(mcvip, genome)

# Saving as RData for fast loading later
save(control_gc, mcvip_gc, file = paste0("data/GpC/", tissue, "/GpC_new.RData"))

# Loading if needed
load(paste0("data/GpC/", tissue, "/GpC_new.RData"))

########
# DMRs #
########

# Load DMRs
load(paste0("data/dmrs/", tissue, "/CG.RData"))
load(paste0("data/dmrs/", tissue, "/CHG.RData"))
load(paste0("data/dmrs/", tissue, "/CHH.RData"))

# Store as list
if (tissue == "leaf") {
    dmrs <- list()
    dmrs[[1]] <- unique(leaf_24h_CpG_ME_CE$"0.15"$"700")
    dmrs[[2]] <- unique(leaf_24h_CHG_ME_CE$"0.15"$"700")
    dmrs[[3]] <- unique(leaf_24h_CHH_ME_CE$"0.15"$"700")
    names(dmrs) <- c("CG", "CHG", "CHH")
} else {
    dmrs <- list()
    dmrs[[1]] <- unique(root_12h_CpG_ME_CE$"0.15"$"700")
    dmrs[[2]] <- unique(root_12h_CHG_ME_CE$"0.15"$"700")
    dmrs[[3]] <- unique(root_12h_CHH_ME_CE$"0.15"$"700")
    names(dmrs) <- c("CG", "CHG", "CHH")
}

# Rename seqlevels
for (i in 1:3) {
    seqlevels(dmrs[[i]]) <- str_remove(seqlevels(dmrs[[i]]), "Chr")
}

# Intersecting ranges
if (pipeline == "split") {
    dmrs <- bind_ranges(dmrs) %>%
        filter(context != "CHH")
} else {
    dmrs <- intersect_ranges(dmrs[[1]], dmrs[[2]]) %>%
        mutate(context = "CG")
}

# Converting DMRs to bed file for TF identification
export_dmr_bed <- FALSE
if (export_dmr_bed) {
    dmrs_bed <- dmrs %>%
        as_tibble() %>%
        mutate(seqnames = paste0("Chr", seqnames), dmr_id = paste0(seqnames, "_", start, "_", end, "_", "merged_", seq_len(length(.))), not_used = NA) %>%
        select(seqnames, start, end, dmr_id, not_used, strand)
    # Exporting intersected DMRs for TF identification
    write_delim(dmrs_bed, paste0("data/dmrs/", tissue, "/merged.bed"), col_names = FALSE, delim = "\t")
} else {
    dmrs_bed <- dmrs %>%
        as_tibble() %>%
        mutate(seqnames = paste0("Chr", seqnames), dmr_id = paste0(seqnames, "_", start, "_", end, "_", "merged_", seq_len(length(.))), not_used = NA) %>%
        select(seqnames, start, end, dmr_id, not_used, strand)
}

####################
# iNOME-seq Signal #
####################

# Adding methylation proportion to GpC data
control_gc_prop <- compute_proportion(control_gc, 4)
mcvip_gc_prop <- compute_proportion(mcvip_gc, 4)

# Removing CHH
control_gc_prop <- control_gc_prop %>%
    filter(context != "CHH")
mcvip_gc_prop <- mcvip_gc_prop %>%
    filter(context != "CHH")

# Exporting methylation as bw
export_gc_prop <- FALSE
if (export_gc_prop) {
    export_cx_report(control_gc_prop, paste0("data/GpC/", tissue, "/"), "control")
    export_cx_report(mcvip_gc_prop, paste0("data/GpC/", tissue, "/"), "mcvip")
}

# Calculate iNOME-seq signal
signal <- compute_signal(control_gc_prop, mcvip_gc_prop)

# Exporting iNOME-seq signal
export_iNOME_signal <- FALSE
if (export_iNOME_signal) {
    export_cx_report(signal, paste0("data/signal/", tissue, "/"), "iNOME_signal")
}

# Merging cytosines
if (pipeline != "split") {
    signal <- signal %>%
        mutate(context = "CG")
}

############################
# Checking cutoff for gaps #
############################

# Calculating cufoffs and plotting data
if (pipeline == "split") {
    cutoffs <- plot_cutoff(dmrs, signal,
        path_for_output = paste0("figures/split/", tissue, "/"),
        prefix = "cutoffs"
    )
} else {
    cutoffs <- plot_cutoff(dmrs, signal,
        path_for_output = paste0("figures/merged/", tissue, "/"),
        prefix = "cutoffs"
    )
}

# Removing plot for simplicity
cutoffs <- cutoffs[[1]]

###########################
# Gaps (length histogram) #
###########################

# Extracting gaps between cytosines (split by context)
if (pipeline == "split") {
    directory <- "split"
} else {
    directory <- "merged"
}

# Importing cpg island locations
cpg_islands <- import("reference/cpg_islands/cpg_islands.gff") %>%
    dropSeqlevels(c("ChrC", "ChrM"), "coarse")
seqlevels(cpg_islands) <- str_remove(seqlevels(cpg_islands), "Chr")

split_gaps <- list()
for (i in seq_len(nrow(cutoffs))) {
    # Extracting gaps with chosen cutoffs
    split_gaps[[i]] <- extract_gaps(
        dmrs, signal,
        pull(cutoffs, cutoff)[i], pull(cutoffs, context)[i],
        cpg_islands, 1
    )
    # Plotting gaps as histogram
    plot_gaps_histogram(
        split_gaps[[i]]$gaps, 125, paste0("figures/", directory, "/", tissue, "/"),
        paste0(pull(cutoffs, context)[i], "_gaps_hist")
    )
    # Merging gaps that are close to eachother
    split_gaps[[i]] <- merge_gaps(split_gaps[[i]])
    # Filtered merge
    split_gaps[[i]]$filtered_gaps <- split_gaps[[i]]$filtered_gaps %>%
        filter(widths <= 60, widths >= 10)
    # Exporting gaps to gff3 for visualisation
    export_gaps(
        split_gaps[[i]]$filtered_gaps, paste0("data/gaps/", tissue, "/"),
        paste0(pull(cutoffs, context)[i], "_", pipeline, "_", pipeline_id, "_gaps")
    )
}
names(split_gaps) <- pull(cutoffs, context)

# Extracting filtered gaps for future analysis
filtered_gaps <- list()
for (i in seq_len(length(split_gaps))) {
    filtered_gaps[[i]] <- split_gaps[[i]]$filtered_gaps
}
names(filtered_gaps) <- names(split_gaps)

##############################
# Gaps (Annotating to genes) #
##############################

# Splitting gff into different features
split_features <- split_gff(
    "reference/TAIR10_GFF3_genes_transposons.gff",
    c("ChrC", "ChrM"), 1:5
)

# Annotating gaps to genes
annotated_gaps_genes <- list()
for (i in seq_len(length(split_gaps))) {
    annotated_gaps_genes[[i]] <- annotate_gaps(
        filtered_gaps[[i]],
        split_features$genes
    )
}
names(annotated_gaps_genes) <- names(split_gaps)

# Saving annotated gaps
save(annotated_gaps_genes, file = paste0(
    "data/gaps/", tissue,
    "/annotated_gaps_", pipeline_id, ".RData"
))

# Filter gaps based on distance to gene
gaps_near_genes <- list()
for (i in seq_len(length(split_gaps))) {
    gaps_near_genes[[i]] <- annotated_gaps_genes[[i]] %>%
        filter(distanceToSite <= 500) %>%
        as_tibble() %>%
        pull(gap)
    gaps_near_genes[[i]] <- filtered_gaps[[i]] %>%
        filter(name %in% gaps_near_genes[[i]])
}
names(gaps_near_genes) <- names(split_gaps)

#######################################
# Gaps (Annotating to split features) #
#######################################

# Annotating gaps to split features
annotated_gaps_features <- list()
annotated_dmrs_features <- list()
for (i in seq_len(length(filtered_gaps))) {
    annotated_gaps_features[[i]] <- annotate_gaps_to_features(
        filtered_gaps[[i]],
        split_features,
        genome
    ) %>%
        mutate(context = names(filtered_gaps)[i])
    annotated_dmrs_features[[i]] <- annotate_gaps_to_features(
        filter(dmrs, context == names(filtered_gaps)[i]),
        split_features,
        genome
    ) %>%
        mutate(context = names(filtered_gaps)[i])
}
annotated_gaps_features <- bind_rows(annotated_gaps_features)
annotated_dmrs_features <- bind_rows(annotated_dmrs_features)

# Replacing expected in gaps with dmrs
annotated_gaps_features <- annotated_gaps_features %>%
    filter(condition != "Expected")
annotated_dmrs_features <- annotated_dmrs_features %>%
    filter(condition != "Expected") %>%
    mutate(condition = "Expected")
annotated_gaps_features <- annotated_gaps_features %>%
    bind_rows(annotated_dmrs_features)

# Plotting annotated gaps to split features
colour_palette2 <- c(
    "#ff3200", "#e9a17c", "#e9e4a6", "#1bb6af", "#0076bb", "#172869"
)
plot_annotated_gaps(
    annotated_gaps_features, paste0("figures/", directory, "/", tissue, "/"),
    "annotated_gaps", colour_palette2
)

# Preparing heatmap comparing TF the expected
tfs <- annotated_gaps_features %>%
    filter(condition == "TF") %>%
    rename(observed = Percentage)
expected <- annotated_gaps_features %>%
    filter(condition == "Expected") %>%
    select(Feature, context, Percentage) %>%
    rename(expected = Percentage)
heatmap_data <- tfs %>%
    inner_join(expected, by = c("Feature" = "Feature", "context" = "context")) %>%
    mutate(score = log(observed / expected)) %>%
    select(Feature, score, context) %>%
    pivot_wider(names_from = context, values_from = score)
feature_ids <- heatmap_data$Feature

heatmap_data <- as.matrix(heatmap_data[, -1])
rownames(heatmap_data) <- feature_ids

if (pipeline == "merged") {
    colnames(heatmap_data) <- "Merged"
}

# Plotting heatmap
heatmap_features <- Heatmap(heatmap_data,
    name = "Log(Observed/Expected)", col = colorRamp2(c(-2, 0, 0.5), c("blue", "white", "red")),
    cluster_rows = FALSE, cluster_columns = FALSE, show_row_names = TRUE,
    show_row_dend = FALSE, show_column_dend = FALSE,
    # column_title = paste0("Log(Observed/Expected)"),
    column_names_gp = grid::gpar(fontsize = 8),
    border_gp = gpar(col = "black")
)
heatmap_features <- draw(heatmap_features)

# Saving heatmaps
pdf(paste0("figures/", directory, "/", tissue, "/", "log_observed_expected.pdf"), width = 10, height = 10)
heatmap_features
dev.off()

###########################################
# Gaps (Distribution of TFs around genes) #
###########################################

# Computing distribution of gaps around TSS
normalised_gaps_tss <- list()
normalised_signal_tss <- list()
for (i in seq_len(length(split_gaps))) {
    # Extracting TSS
    tss <- split_features$genes %>%
        anchor_5p() %>%
        mutate(width = 1)
    # Normalising gaps around TSS
    normalised_gaps_tss[[i]] <- normalise_to_matrix(
        gaps_near_genes[[i]], tss,
        extend = c(500, 1000), bed = TRUE
    )
    # Normalising signal around TSS
    normalised_signal_tss[[i]] <- normalise_to_matrix(
        filter(signal, context == names(gaps_near_genes)[i]), tss,
        extend = c(500, 1000), bed = TRUE
    )
    # Plotting distribution of gaps around TSS
    plot_distribution_tss(normalised_gaps_tss[[i]],
        normalised_signal_tss[[i]],
        extend = c(500, 1000),
        path_for_output = paste0("figures/", directory, "/", tissue, "/"),
        prefix = paste0(names(split_gaps)[i], "_tss_gaps_distribution")
    )
}

# Computing distribution of gaps around genes
normalised_gaps_genes <- list()
normalised_signal_genes <- list()
for (i in seq_len(length(split_gaps))) {
    # Normalising gaps around gene
    normalised_gaps_genes[[i]] <- normalise_to_matrix(
        signal = gaps_near_genes[[i]],
        target = split_features$genes,
        extend = c(500, 500), bed = TRUE
    )
    # Normalising gaps around gene
    normalised_signal_genes[[i]] <- normalise_to_matrix(
        signal = filter(signal, context == names(gaps_near_genes)[i]),
        target = split_features$genes,
        extend = c(500, 500), bed = TRUE
    )
    # Plotting distribution of gaps around genes
    plot_distribution_region(normalised_gaps_genes[[i]],
        normalised_signal_genes[[i]],
        extend = c(500, 500),
        path_for_output = paste0("figures/", directory, "/", tissue, "/"),
        prefix = paste0(names(split_gaps)[i], "_genes_gaps_distribution")
    )
}

#################
# Overlap DNase #
#################

# Loading DNase peaks
dnase_peaks <- read_narrow_peak(
    "data/dnase/DNase_q005_peaks.narrowPeak",
    FALSE
) %>%
    unique() %>%
    keepSeqlevels(as.character(1:5), "coarse")

#  Overlapping gaps with DNase peaks then plotting as a bar chart
overlapped_data <- overlap_data(
    filtered_gaps, dnase_peaks, "Other", paste0("figures/", directory, "/", tissue, "/"),
    "dnase_overlap"
)

#  Overlapping gaps near genes with DNase peaks then plotting as a bar chart
overlapped_data <- overlap_data(
    gaps_near_genes, dnase_peaks, "Other", paste0("figures/", directory, "/", tissue, "/"),
    "dnase_overlap_near_genes"
)

###############################
# Gaps (DNase TF correlation) #
###############################

# Computing correlation between genes and TFs
dnase_TF_correlation <- list()
for (i in seq_len(length(split_gaps))) {
    dnase_TF_correlation[[i]] <-
        compute_correlation(
            gaps_near_genes[[i]],
            dnase_peaks
        ) %>%
        mutate(context = names(split_gaps)[i])
}
dnase_TF_correlation <- bind_rows(dnase_TF_correlation)

if (pipeline == "merged") {
    dnase_TF_correlation <- dnase_TF_correlation %>%
    mutate(context = "Merged")
}

# Plotting gene correlation
plot_gene_gap_correlation(dnase_TF_correlation,
    path_for_output = paste0("figures/", directory, "/", tissue, "/"),
    prefix = "correlation", x_axis_title = "% of DNase peaks per "
)

#################
# Overlap PolII #
#################

# Loading DNase peaks
polII_peaks <- read_narrow_peak(
    "data/polII/SRR5313792_q005_peaks.narrowpeak",
    FALSE
) %>%
    unique() %>%
    keepSeqlevels(as.character(1:5), "coarse")

#  Overlapping gaps with DNase peaks then plotting as a bar chart
overlapped_data_pol <- overlap_data(
    filtered_gaps, polII_peaks, "Other", paste0("figures/", directory, "/", tissue, "/"),
    "polII_overlap"
)

#  Overlapping gaps near genes with DNase peaks then plotting as a bar chart
overlapped_data_pol <- overlap_data(
    gaps_near_genes, polII_peaks, "Other", paste0("figures/", directory, "/", tissue, "/"),
    "polII_overlap_near_genes"
)

###############################
# Gaps (PolII TF correlation) #
###############################

# Computing correlation between genes and TFs
polII_TF_correlation <- list()
for (i in seq_len(length(split_gaps))) {
    polII_TF_correlation[[i]] <-
        compute_correlation(
            filtered_gaps[[i]],
            polII_peaks
        ) %>%
        mutate(context = names(split_gaps)[i])
}
polII_TF_correlation <- bind_rows(polII_TF_correlation)

if (pipeline == "merged") {
    polII_TF_correlation <- polII_TF_correlation %>%
    mutate(context = "Merged")
}

# Plotting gene correlation
plot_gene_gap_correlation(polII_TF_correlation,
    path_for_output = paste0("figures/", directory, "/", tissue, "/"),
    prefix = "correlation_polII",
    x_axis_title = "% of PolII peaks per "
)

#################################################
# Annotating genes with potential PolII pausing #
#################################################

# Annotating gaps to genes
annotated_gaps_polII <- list()
for (i in seq_len(length(split_gaps))) {
    annotated_gaps_polII[[i]] <- annotate_gaps(
        filtered_gaps[[i]],
        split_features$genes,
        bindingRegion = c(0, 200)
    )
}
names(annotated_gaps_polII) <- names(split_gaps)

# Saving annotated gaps
save(annotated_gaps_polII, file = paste0(
    "data/gaps/", tissue,
    "/annotated_gaps_polII_", pipeline_id, ".RData"
))

# Writing gene ids to txt file for GO analysis
for (i in seq_len(length(annotated_gaps_polII))) {
    annotated_gaps_polII[[i]] %>%
        as_tibble() %>%
        select(gene) %>%
        write_delim(paste0(
            "data/gaps/", tissue,
            "/potential_polII_ids_", names(annotated_gaps_polII)[i], "_", pipeline_id, ".txt"
        ), delim = "\t", col_names = FALSE)
}

######################################################
# Plotting go analysis of genes with gaps downstream #
######################################################

if (pipeline == "split") {
    # Loading go data
polII_cg_go <- read_delim(paste0("data/gaps/", tissue, "/potential_polII_ids_CG_S_GO.txt")) %>%
    mutate(query = (queryitem / querytotal) * 100, background = (bgitem / bgtotal) * 100) %>%
    arrange(pvalue) %>%
    slice_head(n = 20) %>%
    arrange(query) %>%
    select(GO_acc, term_type, Term, query, background) %>%
    pivot_longer(-1:-3, names_to = "category", values_to = "percentage")

polII_chg_go <- read_delim(paste0("data/gaps/", tissue, "/potential_polII_ids_CHG_S_GO.txt")) %>%
    mutate(query = (queryitem / querytotal) * 100, background = (bgitem / bgtotal) * 100) %>%
    arrange(pvalue) %>%
    slice_head(n = 20) %>%
    arrange(query) %>%
    select(GO_acc, term_type, Term, query, background) %>%
    pivot_longer(-1:-3, names_to = "category", values_to = "percentage")

# Reordering terms
polII_cg_go$Term <- factor(polII_cg_go$Term, levels = unique(c(polII_cg_go$Term)))
polII_chg_go$Term <- factor(polII_chg_go$Term, levels = unique(c(polII_chg_go$Term)))

# Ploting go analysis
polII_cg_go_plot <- ggplot(polII_cg_go, aes(y = Term, x = percentage, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("#000000", "#993333")) +
    labs(x = "Abundance (%)", y = NULL, fill = NULL, title = "CG") +
    theme_cowplot() +
    theme(axis.line = element_blank(), panel.border = element_rect(colour = "black"))

polII_chg_go_plot <- ggplot(polII_chg_go, aes(y = Term, x = percentage, fill = category)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_manual(values = c("#000000", "#993333")) +
    labs(x = "Abundance (%)", y = NULL, fill = NULL, title = "CHG") +
    theme_cowplot() +
    theme(axis.line = element_blank(), panel.border = element_rect(colour = "black"))

# Combining plots
combined_plots <- polII_cg_go_plot + polII_chg_go_plot + plot_layout(guides = "collect")

# Saving plots
ggsave(paste0("figures/", directory, "/", tissue, "/polII_geens_go_analysis.pdf"), combined_plots, height = 10, width = 10)
}

#################################
# Preparing gaps for meme suite #
#################################

# Selecting dmrs that do not have any gaps inside
dmrs_no_gaps <- vector("list", length(filtered_gaps))
names(dmrs_no_gaps) <- names(filtered_gaps)
for (i in seq_len(length(filtered_gaps))) {
    dmrs_no_gaps[[i]] <- subsetByOverlaps(
        filter(dmrs, context == names(filtered_gaps)[i]), filtered_gaps[[i]],
        invert = TRUE
    )
    # Reducing the ranges of the dmrs without gaps
    dmrs_no_gaps[[i]] <- reduce_ranges(dmrs_no_gaps[[i]])
    dmrs_no_gaps[[i]] <- dmrs_no_gaps[[i]][ width(dmrs_no_gaps[[i]]) > 10]
    # Changing seq levels to match fasta file
    seqlevels(dmrs_no_gaps[[i]]) <- paste0("Chr", seqlevels(dmrs_no_gaps[[i]]))
    # Exporting dmrs without gaps to bed file
    export(dmrs_no_gaps[[i]], paste0(
        "data/gaps/", tissue,
        "/dmrs_no_gaps_", names(dmrs_no_gaps)[i], "_", pipeline_id, ".bed"
    ))
}

# Exporting filtered gaps as bed file
for (i in seq_len(length(filtered_gaps))) {
    # Changing seq levels to match fasta file
    seqlevels(filtered_gaps[[i]]) <- paste0("Chr", seqlevels(filtered_gaps[[i]]))
    export(filtered_gaps[[i]], paste0(
        "data/gaps/", tissue,
        "/filtered_gaps_", names(dmrs_no_gaps)[i], "_", pipeline_id, ".bed"
    ))
}

# use this to get fasta files for gaps and background
# for i in *bed; do bedtools getfasta -fi TAIR10_chr_all.fasta -bed $i -fo ${i/bed/fasta}; done

#######################
# Quantifying Overlap #
#######################

# Storing Chromosome lengths
chromosome_lengths <- data.frame(chr = names(genome), length = width(genome)) %>%
    mutate(seqnames = chr, start = 1, end = length) %>%
    select(seqnames, start, end) %>%
    makeGRangesFromDataFrame(keep.extra.columns = TRUE)

# Plotting overlap with DNase
dnase_test <- vector("list", length(gaps_near_genes))
names(dnase_test) <- names(gaps_near_genes)
for (i in seq_len(length(gaps_near_genes))) {
    dnase_test[[i]] <- permTest(
        A = gaps_near_genes[[i]], B = dnase_peaks,
        randomize.function = randomizeRegions,
        evaluate.function = numOverlaps,
        genome = chromosome_lengths,
        ntimes = 1000
    )
    dnase_perm_plot <- as.ggplot(function() plot(dnase_test[[i]]))
    ggsave(plot = dnase_perm_plot, filename = paste0(
        "figures/",
        directory, "/", tissue, "/dnase_perm_test_", names(dnase_test)[i], ".pdf"
    ), height = 10, width = 10)
}

# Plotting overlap with polII
polII_test <- vector("list", length(gaps_near_genes))
names(polII_test) <- names(gaps_near_genes)
for (i in seq_len(length(gaps_near_genes))) {
    polII_test[[i]] <- permTest(
        A = gaps_near_genes[[i]], B = polII_peaks,
        randomize.function = randomizeRegions,
        evaluate.function = numOverlaps,
        genome = chromosome_lengths,
        ntimes = 1000
    )
    polII_perm_plot <- as.ggplot(function() plot(polII_test[[i]]))
    ggsave(plot = polII_perm_plot, filename = paste0(
        "figures/",
        directory, "/", tissue, "/polII_perm_test_", names(polII_test)[i], ".pdf"
    ), height = 10, width = 10)
}

###########
# RNA-seq #
###########

# Calculating gene sizes
calculate_gene_sizes <- FALSE
if (calculate_gene_sizes) {
    GFF <- makeTxDbFromGFF("reference/TAIR10_GFF3_genes_transposons.gff")
    seqlevels(GFF) <- c(1:5, "C", "M")
    GFF <- keepSeqlevels(GFF, c(1:5), "coarse")
    exons_list_per_gene <- exonsBy(GFF, by = "gene")
    exonic_gene_sizes <- list()
    for (gene in names(exons_list_per_gene)) {
        exonic_gene_sizes[[gene]] <- sum(width(
            reduce_ranges(exons_list_per_gene[[gene]])
        ))
    }
    lengths <- tibble(
        id = names(exonic_gene_sizes),
        length = unlist(exonic_gene_sizes), row.names = NULL
    )
    save(lengths, file = "data/expression/4_stats/gene_lengths.RData")
}

# Load Gene lengths
load("data/expression/4_stats/gene_lengths.RData")

# Load RNA-seq counts, calculate tpm then catagorise genes based on expression
if (tissue == "leaf") {
    leaf_expression <- read_delim(
        "data/expression/4_stats/leaf_ReadsPerGene.out.tab",
        col_names = FALSE
    ) %>%
        select(1:2) %>%
        rename(id = X1, count = X2) %>%
        inner_join(lengths)
    gene_id <- pull(leaf_expression, id)

    leaf_expression <- tpm(leaf_expression[-1], leaf_expression$length) %>%
        as_tibble() %>%
        select(-length, tpm = count) %>%
        mutate(id = gene_id) %>%
        relocate(id)
    # Saving TPMs
    write_delim(leaf_expression, "data/expression/4_stats/leaf_tpm.csv",
        delim = ","
    )
} else {
    root_expression <- read_delim(
        "data/expression/4_stats/root_ReadsPerGene.out.tab",
        col_names = FALSE
    ) %>%
        select(1:2) %>%
        rename(id = X1, count = X2) %>%
        inner_join(lengths)
    gene_id <- pull(root_expression, id)
    root_expression <- tpm(root_expression[-1], root_expression$length) %>%
        as_tibble() %>%
        select(-length, tpm = count) %>%
        mutate(id = gene_id) %>%
        relocate(id)
    # Saving TPMs
    write_delim(root_expression, "data/expression/4_stats/root_tpm.csv",
        delim = ","
    )
}

#####################################
# Annotating genes with CpG islands #
#####################################

# Importing cpg island locations
cpg_islands <- import("reference/cpg_islands/cpg_islands.gff") %>%
    dropSeqlevels(c("ChrC", "ChrM"), "coarse") %>%
    select(name = ID)
seqlevels(cpg_islands) <- str_remove(seqlevels(cpg_islands), "Chr")

# Annotating CpG islands to genes
annotated_islands_genes <- annotate_gaps(
    cpg_islands,
    split_features$genes,
    bindingRegion = c(-1000, 1)
)

# Extracting gene ids
annotated_islands_ids <- annotated_islands_genes %>%
    as_tibble() %>%
    pull(gene) %>%
    unique()

# Reannotating gaps closer to genes
annotated_gaps_promoter <- list()
for (i in seq_len(length(filtered_gaps))) {
    temp <- filtered_gaps[[i]]
    seqlevels(temp) <- str_remove(seqlevels(temp), "Chr")
    annotated_gaps_promoter[[i]] <- annotate_gaps(
        temp,
        split_features$genes,
        bindingRegion = c(-1000, 1)
    )
}
names(annotated_gaps_promoter) <- names(filtered_gaps)

#############################################
# Comparing Genes with gap to genes without #
#############################################

# Loading Annotated gaps and gene expression
if (tissue == "leaf") {
    load(paste0("data/gaps/leaf/annotated_gaps_", pipeline_id, ".RData"))
    gene_expression <- read_delim("data/expression/4_stats/leaf_tpm.csv",
        delim = ","
    )
} else {
    load(paste0("data/gaps/root/annotated_gaps_", pipeline_id, ".RData"))
    gene_expression <- read_delim("data/expression/4_stats/root_tpm.csv",
        delim = ","
    )
}

# Categorising genes
tf_expression <- list()
for (i in seq_len(length(annotated_gaps_promoter))) {
    gap_id <- annotated_gaps_promoter[[i]] %>%
        as_tibble() %>%
        pull(gene) %>%
        unique()
    tf_expression[[i]] <- gene_expression %>%
        mutate(tpm = log(tpm + 1)) %>%
        mutate(island = ifelse(id %in% annotated_islands_ids, "yes", "no")) %>%
        mutate(gap = ifelse(id %in% gap_id, "yes", "no")) %>%
        mutate(label = ifelse(island == "no" & gap == "no", "No TF No Island",
            ifelse(island == "yes" & gap == "yes", "TF Island",
                ifelse(island == "no" & gap == "yes", "TF No Island",
                    ifelse(island == "yes" & gap == "no", "No TF Island", NA)
                )
            )
        ))
}
names(tf_expression) <- names(annotated_gaps_promoter)

# Plotting Expression of annotated genes vs other genes
my_comparisons <- list(
    c("No TF Island", "No TF No Island"),
    c("No TF Island", "TF Island"),
    c("No TF Island", "TF No Island"),
    c("No TF No Island", "TF Island"),
    c("No TF No Island", "TF No Island"),
    c("TF Island", "TF No Island")
)

if (pipeline == "merged") {
    names(tf_expression) <-  "Merged"
}

tf_expression_plot <- list()
tf_expression_plot[[1]] <- ggplot(tf_expression[[1]], aes(x = label, y = tpm, fill = label)) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
    labs(x = NULL, y = "log(tpm + 1)", title = names(tf_expression)[1], fill = NULL) +
    theme_cowplot() +
    theme(
        axis.line = element_blank(), panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    )

tf_expression_plot[[2]] <- ggplot(tf_expression[[2]], aes(x = label, y = tpm, fill = label)) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
    labs(x = NULL, y = "log(tpm + 1)", title = names(tf_expression)[2], fill = NULL) +
    theme_cowplot() +
    theme(
        axis.line = element_blank(), panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    )

tf_expression_plot[[3]] <- ggplot(tf_expression[[3]], aes(x = label, y = tpm, fill = label)) +
    geom_boxplot(outlier.shape = NA) +
    stat_compare_means(comparisons = my_comparisons, method = "t.test", label = "p.signif") +
    labs(x = NULL, y = "log(tpm + 1)", title = names(tf_expression)[3], fill = NULL) +
    theme_cowplot() +
    theme(
        axis.line = element_blank(), panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    )

if (pipeline == "split") {
    combined_plots <- tf_expression_plot[[1]] + tf_expression_plot[[2]] + tf_expression_plot[[3]] + plot_layout(guides = "collect")
} else {
    combined_plots <- tf_expression_plot[[1]]
}

# Saving plot
ggsave(paste0(
    "figures/",
    directory, "/", tissue, "/gene_expression.pdf"
), combined_plots, width = 15.9, height = 11.7)

# Count genes in each category
tf_expression_count_plot <- list()
tf_expression_count_plot[[1]] <- tf_expression[[1]] %>%
    group_by(label) %>%
    summarise(count = dplyr::n()) %>%
    ggplot(aes(x = label, y = count, fill = label)) +
    geom_bar(stat = "identity", colour = "black") +
    geom_text(aes(label = count), vjust = -0.9, color = "black") +
    labs(x = NULL, y = "Gene Count", title = names(tf_expression)[1], fill = NULL) +
    theme_cowplot() +
    theme(
        axis.line = element_blank(), panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    )

tf_expression_count_plot[[2]] <- tf_expression[[2]] %>%
    group_by(label) %>%
    summarise(count = dplyr::n()) %>%
    ggplot(aes(x = label, y = count, fill = label)) +
    geom_bar(stat = "identity", colour = "black") +
    geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.9, color = "black") +
    labs(x = NULL, y = "Gene Count", title = names(tf_expression)[2], fill = NULL) +
    theme_cowplot() +
    theme(axis.line = element_blank(), panel.border = element_rect(colour = "black"),
    axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5))

tf_expression_count_plot[[3]] <- tf_expression[[3]] %>%
    group_by(label) %>%
    summarise(count = dplyr::n()) %>%
    ggplot(aes(x = label, y = count, fill = label)) +
    geom_bar(stat = "identity", colour = "black") +
    geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.9, color = "black") +
    labs(x = NULL, y = "Gene Count", title = names(tf_expression)[3], fill = NULL) +
    theme_cowplot() +
    theme(
        axis.line = element_blank(), panel.border = element_rect(colour = "black"),
        axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5)
    )

if (pipeline == "split") {
    combined_plots <- tf_expression_count_plot[[1]] + tf_expression_count_plot[[2]] +
        tf_expression_count_plot[[3]] + plot_layout(guides = "collect")
} else {
    combined_plots <- tf_expression_count_plot[[1]]
}

# Saving plot
ggsave(paste0(
    "figures/",
    directory, "/", tissue, "/gene_expression_count.pdf"
), combined_plots, width = 15.9, height = 11.7)

#############################
# Identified TFs expression #
#############################

# Loading Annotated gaps and gene expression
if (tissue == "leaf") {
    load(paste0("data/gaps/leaf/annotated_gaps_", pipeline_id, ".RData"))
    gene_expression <- read_delim("data/expression/4_stats/leaf_tpm.csv",
        delim = ","
    )
} else {
    load(paste0("data/gaps/root/annotated_gaps_", pipeline_id, ".RData"))
    gene_expression <- read_delim("data/expression/4_stats/root_tpm.csv",
        delim = ","
    )
}

# Importing identified TFs and exporting ids that need conversion
if (pipeline == "merged") {
    contexts <- "merged"
} else {
    contexts <- c("CG", "CHG", "CHH")
}

tf_ids <- list()
for (i in contexts) {
    tf_ids[[i]] <- read_delim(paste0("D:/iNOME/data/gaps/", tissue, "/", i, "/xstreme.tsv"), delim = "\t") %>%
        select(ALT_ID) %>%
        distinct() %>%
        filter(str_detect(ALT_ID, "STREME", negate = TRUE))
    temp <- tf_ids[[i]] %>%
        filter(str_detect(ALT_ID, "[aA][tT][0-9][gG]", negate = TRUE))
    write_delim(temp,
        paste0("D:/iNOME/data/gaps/", tissue, "/", i, "/OFFICIAL_GENE_SYMBOL.tsv"),
        col_names = FALSE
    )
}

# Importing corrected gene ids then exporting ids that failed conversion to 
# convert manually
for (i in contexts) {
    file_to_import <- dir(paste0("D:/iNOME/data/gaps/",tissue,"/", i), pattern = "gProfiler", full.names = TRUE)
    corrected_ids <- read_delim(file_to_import) %>%
        select(initial_alias, converted_alias)
    id_that_need_correcting <- corrected_ids %>%
        filter(converted_alias == "None")
    write_delim(id_that_need_correcting, paste0("D:/iNOME/data/gaps/", tissue, "/", i, "/manual_correction.txt"), col_names = FALSE)
}

# Combining automatic corrected ids and manually corrected ids
all_ids <- list()
for (i in contexts) {
    file_to_import <- dir(paste0("D:/iNOME/data/gaps/", tissue, "/", i), pattern = "gProfiler", full.names = TRUE)
    file_to_import2 <- dir(paste0("D:/iNOME/data/gaps/", tissue, "/", i), pattern = "corrected", full.names = TRUE)
    auto_corrected_ids <- read_delim(file_to_import) %>%
        select(initial_alias, converted_alias)
    manual_corrected_ids <- read_delim(file_to_import2, col_names = c("initial_alias", "converted_alias"))
    converted_ids <- auto_corrected_ids %>%
        full_join(manual_corrected_ids, by = c("initial_alias" = "initial_alias")) %>%
        distinct() %>%
        mutate(temp = ifelse(converted_alias.x == "None" & !is.na(converted_alias.y), converted_alias.y,
            ifelse(converted_alias.x != "None", converted_alias.x, "NA")
        )) %>%
        mutate(initial_alias, converted_alias = temp) %>%
            select(initial_alias, converted_alias) %>%
            filter(str_detect(converted_alias, "AT"))
    already_correct_ids <- tf_ids[[i]] %>%
        filter(str_detect(ALT_ID, "[aA][tT][0-9][gG]")) %>%
        select(initial_alias = ALT_ID, converted_alias = ALT_ID)
    all_ids[[i]] <- bind_rows(already_correct_ids, converted_ids) %>%
        mutate(converted_alias = str_to_upper(converted_alias))
}

# Exporting genes that had a correct id but need a gene name
for (i in contexts) {
    all_ids[[i]] %>%
        mutate(initial_alias = str_to_upper(initial_alias)) %>%
        filter(str_detect(initial_alias, "AT[0-9]")) %>%
        select(initial_alias) %>%
        write_delim(paste0("D:/iNOME/data/gaps/", tissue, "/", i, "/tair_gene_ids_to_convert.txt"), delim = "\t", col_names = FALSE)
}

# Importing correct gene name for these ids
all_ids_gene_names <- list()
for (i in contexts) {
    ATG_to_gene_ids <- read_delim(paste0("D:/iNOME/data/gaps/", tissue, "/", i, "/tair_gene_ids_converted.txt"),
        delim = "\t"
    ) %>%
    filter(converted_alias != "None") %>%
    mutate(x1 = initial_alias, x2 = converted_alias, initial_alias = x2, converted_alias = x1) %>%
    select(initial_alias, converted_alias)
    all_ids_gene_names[[i]] <- all_ids[[i]] %>%
        mutate(initial_alias = str_to_upper(initial_alias)) %>%
        filter(str_detect(initial_alias, "AT[0-9]", negate = TRUE)) %>%
        bind_rows(ATG_to_gene_ids)
}

# Creating matrix of expression for identified TFs
rna_counts_heatmap <- vector("list", length(all_ids))
names(rna_counts_heatmap) <- names(all_ids)
for (i in seq_len(length(all_ids))) {
    rna_counts_heatmap[[i]] <- gene_expression %>%
        filter(id %in% all_ids_gene_names[[i]]$converted_alias) %>%
        rename(genes = id) %>%
        filter(tpm > 15)
    gene_ids <- pull(rna_counts_heatmap[[i]], genes)
    tf_ids <- all_ids_gene_names[[i]] %>%
        filter(converted_alias %in% gene_ids) %>%
        mutate(initial_alias = str_to_upper(initial_alias)) %>%
        pull(initial_alias)
    sample_ids <- colnames(rna_counts_heatmap[[i]][-1])
    rna_counts_heatmap[[i]] <- rna_counts_heatmap[[i]] %>%
        select(-genes) %>%
        as.matrix() # %>%
    # t() %>%
    # scale() %>%
    # t()
    colnames(rna_counts_heatmap[[i]]) <- sample_ids
    #rownames(rna_counts_heatmap[[i]]) <- gene_ids
    rownames(rna_counts_heatmap[[i]]) <- tf_ids
}

# Ploitting heatmaps
if (pipeline == "split") {
    CG_heatmap <- Heatmap(rna_counts_heatmap[[1]],
        name = "Expression", col = colorRampPalette(c("white", "red"))(100),
        cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = TRUE,
        show_row_dend = FALSE, row_names_side = "left", row_title = "Genes",
        column_title = paste0(names(rna_counts_heatmap)[1], " TF Expression"),
        row_names_gp = gpar(fontsize = 6)
    )
    CG_heatmap <- draw(CG_heatmap)

    CHG_heatmap <- Heatmap(rna_counts_heatmap[[2]],
        name = "Expression", col = colorRampPalette(c("white", "red"))(100),
        cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = TRUE,
        show_row_dend = FALSE, row_names_side = "left", row_title = "Genes",
        column_title = paste0(names(rna_counts_heatmap)[2], " TF Expression"),
        row_names_gp = gpar(fontsize = 6)
    )
    CHG_heatmap <- draw(CHG_heatmap)

    CHH_heatmap <- Heatmap(rna_counts_heatmap[[3]],
        name = "Expression", col = colorRampPalette(c("white", "red"))(100),
        cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = TRUE,
        show_row_dend = FALSE, row_names_side = "left", row_title = "Genes",
        column_title = paste0(names(rna_counts_heatmap)[3], " TF Expression"),
        row_names_gp = gpar(fontsize = 6)
    )
    CHH_heatmap <- draw(CHH_heatmap)

    # Saving heatmaps
    pdf(paste0(
        "figures/",
        pipeline, "/", tissue, "/CG_TF_expression.pdf"
    ), height = 8, width = 2)
    print(CG_heatmap)
    dev.off()

    pdf(paste0(
        "figures/",
        pipeline, "/", tissue, "/CHG_TF_expression.pdf"
    ), height = 8, width = 2)
    print(CHG_heatmap)
    dev.off()

    pdf(paste0(
        "figures/",
        pipeline, "/", tissue, "/CHH_TF_expression.pdf"
    ), height = 8, width = 2)
    print(CHH_heatmap)
    dev.off()
} else {
    merged_heatmap <- Heatmap(rna_counts_heatmap[[1]],
        name = "Expression", col = colorRampPalette(c("white", "red"))(100),
        cluster_rows = TRUE, cluster_columns = FALSE, show_row_names = TRUE,
        show_row_dend = FALSE, row_names_side = "left", row_title = "Genes",
        column_title = paste0(names(rna_counts_heatmap)[1], " TF Expression"),
        row_names_gp = gpar(fontsize = 6)
    )
    merged_heatmap <- draw(merged_heatmap)
    pdf(paste0(
        "figures/",
        pipeline, "/", tissue, "/merged_TF_expression.pdf"
    ), height = 8, width = 2)
    merged_heatmap
    dev.off()
}

###########################
# Comparing root and leaf #
###########################

# Splitting gff into different features
split_features <- split_gff(
    "reference/TAIR10_GFF3_genes_transposons.gff",
    c("ChrC", "ChrM"), 1:5
)

# Importing Gene expression
leaf_expression <- read_delim("data/expression/4_stats/leaf_tpm.csv",
        delim = ","
    )
root_expression <- read_delim("data/expression/4_stats/root_tpm.csv",
        delim = ","
    )

# Loading Annotated gaps and gene expression
load(paste0("data/gaps/leaf/annotated_gaps_", pipeline_id, ".RData"))
leaf_annotated_gaps <- annotated_gaps_genes
load(paste0("data/gaps/root/annotated_gaps_", pipeline_id, ".RData"))
root_annotated_gaps <- annotated_gaps_genes

# Counting genes in each category and plotting counts
if (pipeline == "merged") {
    contexts <- "merged"
} else {
    contexts <- c("CG", "CHG")
}

categorised_genes_count_plot <- list()
categorised_genes_expression_plot <- list()
for (i in seq_len(length(contexts))) {
    categorised_genes <- split_features$genes %>%
        mutate(category = "Neither") %>%
        mutate(category = ifelse(gene_id %in% leaf_annotated_gaps[[i]]$gene & gene_id %in% root_annotated_gaps[[i]]$gene, "Both", category)) %>%
        mutate(category = ifelse(gene_id %in% leaf_annotated_gaps[[i]]$gene & !gene_id %in% root_annotated_gaps[[i]]$gene, "Leaf Only", category)) %>%
        mutate(category = ifelse(!gene_id %in% leaf_annotated_gaps[[i]]$gene & gene_id %in% root_annotated_gaps[[i]]$gene, "Root Only", category)) %>%
        mutate(category = as.factor(category)) %>%
        as_tibble() %>%
        inner_join(root_expression, by = c("gene_id" = "id")) %>%
        inner_join(leaf_expression, by = c("gene_id" = "id")) %>%
        rename(root = tpm.x, leaf = tpm.y)
    categorised_genes_counts <- categorised_genes %>%
        group_by(category) %>%
        summarise(count = n())
    categorised_genes_count_plot[[i]] <- ggplot(categorised_genes_counts, aes(x = category, y = count, fill = category)) +
        geom_bar(stat = "identity", colour = "black") +
        geom_text(aes(label = count), position = position_dodge(width = 0.9), vjust = -0.9, color = "black") +
        labs(x = "Category", y = "Gene count", fill = NULL) +
        theme_cowplot() +
        theme(axis.line = element_blank(), panel.border = element_rect(colour = "black"))
    categorised_genes <- categorised_genes %>%
        pivot_longer(-seqnames:-category, names_to = "tissue", values_to = "expression") %>%
        mutate(expression = log10(expression + 1))
    my_comparisons <- c("leaf", "root")
    categorised_genes_expression_plot[[i]] <- ggplot(categorised_genes, aes(x = category, y = expression, fill = tissue)) +
        geom_boxplot(outlier.shape = NA) +
        stat_compare_means(aes(group = tissue), method = "t.test", label = "p.signif") +
        theme_cowplot() +
        theme(axis.line = element_blank(), panel.border = element_rect(colour = "black"))
}

if (pipeline == "merged") {
    combined_plots_count <- categorised_genes_count_plot[[1]]
    combined_plots_expression <- categorised_genes_expression_plot[[1]]
} else {
    combined_plots_count <- categorised_genes_count_plot[[1]] + categorised_genes_count_plot[[2]] + plot_layout(guides = "collect")
    combined_plots_expression <- categorised_genes_expression_plot[[1]] + categorised_genes_expression_plot[[2]] + plot_layout(guides = "collect")
}

# Saving plot
ggsave(paste0(
    "figures/",
    directory, "/categorised_genes_count_plot.pdf"
), combined_plots_count, width = 15.9, height = 11.7)

ggsave(paste0(
    "figures/",
    directory, "/categorised_genes_expression_plot.pdf"
), combined_plots_expression, width = 15.9, height = 11.7)
