# Import CX_reports
import_report <- function(path_to_file) {
    output <- read_delim(path_to_file,
        col_names = c(
            "seqnames", "start", "strand",
            "readsM", "readsN", "context",
            "trinucleotide_context"
        )
    ) %>%
        mutate(end = start) %>%
        filter(!seqnames %in% c("Pt", "Mt")) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    return(output)
}

# mport fasta genome
import_genome <- function(path_to_fasta) {
    output <- import(path_to_fasta)
    output <- output[str_detect(names(output), "Chr\\d")]
    names(output) <- str_match(names(output), "Chr(\\d)")[, 2]
    return(output)
}

# Extract sequence information for cytosines from CX report
extract_sequence <- function(methylation_data, fasta, type = "GC") {
    chromosome_regions <- GRanges(
        seqnames = names(fasta),
        ranges = IRanges(start = 2, end = width(fasta) - 1)
    )
    cytosine_positions <- methylation_data %>%
        anchor_3p() %>%
        mutate(width = 2) %>%
        subsetByOverlaps(chromosome_regions)
    if (type == "GC") {
        output <- cytosine_positions[as.character(genome[cytosine_positions]) == "GC"] %>%
            anchor_3p() %>%
            mutate(width = 1)
    } else if (type == "other") {
        output <- cytosine_positions[as.character(genome[cytosine_positions]) != "GC"] %>%
            anchor_3p() %>%
            mutate(width = 1)
    } else {
        stop("Use type = 'GC' or 'other")
    }
    return(output)
}

# Function to compute dmrs
compute_dmrs_bins <- function(control, treatment, contexts,
                              binSize,
                              pValueThreshold = 0.05, minCytosinesCount = 1,
                              minProportionDifference,
                              minReadsPerCytosine = 4, cores = 4,
                              minGap = binSize * 2) {
    output <- list()
    for (context in contexts) {
        output[[context]] <- computeDMRs(control, treatment,
            context = context,
            binSize = binSize, pValueThreshold = pValueThreshold,
            minCytosinesCount = minCytosinesCount,
            minProportionDifference = minProportionDifference,
            minReadsPerCytosine = minReadsPerCytosine, cores = cores,
            minGap = minGap
        )
    }
    output <- bind_ranges(output)
    return(output)
}

# Export dmrs as gff3 for jbrowse
export_dmrs <- function(dmrs, path_for_output, prefix) {
    dmrs_all <- dmrs %>%
        mutate(difference = proportion2 - proportion1) %>%
        dplyr::select(type = context, score = difference)
    export(dmrs_all, paste0(path_for_output, prefix, "_all.gff3"))
    for (context in unique(dmrs_all$type)) {
        dmrs_all %>%
            filter(type == context) %>%
            export(paste0(path_for_output, prefix, "_", context, ".gff3"))
    }
}

# Compute methylation proportion
compute_proportion <- function(report, min_reads = 4) {
    output <- report %>%
        mutate(score = (readsM / readsN)) %>%
        filter(readsM + readsN >= min_reads) %>%
        dplyr::select(context, score) %>%
        filter(!is.na(score))
    return(output)
}

# Compute iNOME-seq signal from CX reports
compute_signal <- function(control, treatment) {
    control <- control %>%
        mutate(control_score = score) %>%
        dplyr::select(context, control_score) %>%
        as_tibble()
    treatment <- treatment %>%
        mutate(treatment_score = score) %>%
        dplyr::select(context, treatment_score) %>%
        as_tibble()
    output <- control %>%
        inner_join(treatment, by = c(
            "seqnames" = "seqnames", "start" = "start",
            "strand" = "strand", "context" = "context"
        )) %>%
        mutate(score = treatment_score - control_score) %>%
            filter(!is.na(score)) %>%
        dplyr::select(seqnames, start, end = end.x, strand, context, score) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    return(output)
}

# Export CX report
export_cx_report <- function(report, path_to_file, prefix) {
    for (i in unique(report$context)) {
        report %>%
            filter(context == i) %>%
            "seqlengths<-"(end(getWholeChromosomes(.))) %>%
            dplyr::select(-context) %>%
            export(paste0(path_to_file, prefix, "_", i, ".bw"))
    }
}

# Calculate and plot cutoffs for gap identification
plot_cutoff <- function(dmrs, signal, contexts = c("CG", "CHG", "CHH"), path_for_output = NULL, prefix = NULL) {
    overlaps <- c()
    for (i in contexts) {
        temp_dmrs <- filter(dmrs, context == i)
        temp_signal <- filter(signal, context == i)
        overlaps <- bind_ranges(
            overlaps,
            subsetByOverlaps(temp_signal, temp_dmrs)
        )
    }
    output <- list()
    raw <- overlaps %>%
        filter(score >= 0) %>%
        as_tibble() %>%
        dplyr::select(context, score)
    output[[1]] <- raw %>%
        group_by(context) %>%
        summarise(cutoff = mean(score))
    output[[2]] <- ggplot(raw, aes(x = context, y = score, fill = context)) +
        geom_boxplot() +
        geom_bar(
            data = output[[1]], aes(x = context, y = cutoff, fill = context),
            stat = "identity", alpha = 0.2
        ) +
        labs(
            x = "Methylation Context", y = "Methylation Proportion",
            fill = "Context"
        ) +
        theme_cowplot() +
        grids(linetype = "dashed") +
        theme(
            axis.line = element_blank(),
            panel.border = element_rect(colour = "black")
        )
    if (!is.null(path_for_output) > 0 && !is.null(prefix)) {
        ggsave(
            filename = paste0(path_for_output, prefix, ".pdf"),
            plot = output[[2]],
            bg = "white", height = 10, width = 10
        )
    }
    names(output) <- c("Cutoffs", "Plot")
    return(output)
}

# Extract gaps and filter them based on size
extract_gaps <- function(dmrs, signal, meth_cutoff, filter_context, cpg_islands, cytosines_in_gaps = 1) {
    signal <- filter(signal, context == filter_context, score >= meth_cutoff)
    dmrs <- dmrs %>%
        filter(context == filter_context) %>%
        mutate(
            label = paste0(seqnames, ":", start(.), "-", end(.))
        ) %>%
        dplyr::select(type = context, label)
    output <- list()
    output[[1]] <- subsetByOverlaps(signal, dmrs) %>%
        gaps() %>%
        mutate(widths = width(.), strand = "*") %>%
        filter(widths <= 1000)
    index <- findOverlaps(signal, output[[1]], type = "within") %>%
        as_tibble() %>%
        group_by(subjectHits) %>%
        summarise(count = dplyr::n()) %>%
        filter(count >= cytosines_in_gaps) %>%
        pull(subjectHits)
    output[[1]] <- output[[1]][index]
    if (!is.null(cpg_islands)) {
        index <- findOverlaps(cpg_islands, output[[1]]) %>%
            subjectHits() %>%
            unique()
        output[[1]] <- output[[1]][-index]
    }
    output[[2]] <- filter(output[[1]], widths <= 60, widths >= 10)
    names(output) <- c("gaps", "filtered_gaps")
    return(output)
}

# Merge close gaps
merge_gaps <- function(gaps_output) {
    output <- vector("list", 2)
    names(output) <- c("gaps", "filtered_gaps")
    output[[1]] <- gaps_output[[1]] %>%
        stretch(2) %>%
        reduce_ranges() %>%
        stretch(-2) %>%
        mutate(widths = width(.), name = paste0("gap_", seq(1, length(.))))
    output[[2]] <- gaps_output[[2]] %>%
        stretch(2) %>%
        reduce_ranges() %>%
        stretch(-2) %>%
        mutate(widths = width(.), name = paste0("gap_", seq(1, length(.))))
    return(output)
}

# Plot histogram of gaps
plot_gaps_histogram <- function(gaps, bins = 125, path_for_output = NULL,
                                prefix = NULL) {
    output <- gaps %>%
        as_tibble() %>%
        filter(widths > 1) %>%
        ggplot(aes(x = widths)) +
        geom_histogram(aes(y = ..density..),
            bins = bins, colour = "black",
            fill = "white"
        ) +
        geom_density(colour = 4, fill = 4, alpha = 0.25) +
        geom_vline(
            xintercept = c(50, 150, 300, 450), lty = "dashed",
            colour = "red"
        ) +
        labs(x = "Length distribution of unmethylated runs", y = "Count") +
        xlim(0,1000) +
        theme_cowplot() +
        theme(panel.border = element_rect(colour = "black"))
    if (!is.null(path_for_output) > 0 && !is.null(prefix)) {
        ggsave(
            filename = paste0(path_for_output, prefix, ".pdf"), plot = output,
            bg = "white", height = 10, width = 10
        )
    }
    return(output)
}

# export gaps as gff3
export_gaps <- function(gaps, path_for_output, prefix) {
    gaps %>%
        mutate(type = "gap", score = widths) %>%
        export(paste0(path_for_output, prefix, ".gff3"))
}

# Loading annoation and extracting genes
prepare_genes <- function(path_to_annotation, seq_levels = as.character(1:5)) {
    output <- import(path_to_annotation) %>%
        filter(type == "gene") %>%
        dropSeqlevels(c("ChrC", "ChrM"), "coarse")
    seqlevels(output) <- seq_levels
    return(output)
}

# Annotate gaps to genes
annotate_gaps <- function(gaps, genes,
                          bindingType = "startSite",
                          bindingRegion = c(-3000, 1)) {
    output <- annoPeaks(gaps, genes,
        bindingType = bindingType,
        bindingRegion = bindingRegion
    ) %>%
        dplyr::select(
            gap = name, feature.ranges, gene = gene_id, distanceToSite,
            location = insideFeature
        ) %>%
        arrange(distanceToSite)
    return(output)
}

# Normalise data to matrix then smooth line
slider_median <- function(x) {
    slide_dbl(x, ~ median(.x, na.rm = TRUE), .before = 3, .after = 3)
}
slider_mean <- function(x) {
    slide_dbl(x, ~ mean(.x, na.rm = TRUE), .before = 3, .after = 3)
}

normalise_to_matrix <- function(signal, target, extend = 2000, w = max(extend) / 50,
                                meth = FALSE, bed = FALSE,
                                value_column = if (!bed) {
                                    "score"
                                }, # trick to pass a NULL if bed is TRUE
                                mapping_column = NULL, background = NA,
                                mean_mode = case_when(meth ~ "absolute", bed ~ "coverage", TRUE ~ "w0"),
                                include_target = any(width(target) > 1),
                                target_ratio = median(width(target)) / (sum(extend) + median(width(target)) + sum(extend)),
                                k = min(c(20, min(width(target)))), smooth = FALSE, smooth_fun = slider_median,
                                keep = c(0, 0.99), flip_upstream = FALSE) {
        normalizeToMatrix(
            signal = signal, target = target, extend = extend, w = w, value_column = value_column,
            mapping_column = mapping_column, background = background, mean_mode = mean_mode,
            include_target = include_target, target_ratio = target_ratio, k = k,
            smooth = smooth, smooth_fun = smooth_fun, keep = keep, flip_upstream = flip_upstream
        )
}

# Plotting gaps around genes
plot_distribution_region <- function(normalised_gaps,
                                     normalised_signal,
                                     extend = c(500, 500), label = NULL, 
                                     y_title = NULL, path_for_output = NULL, 
                                     prefix = NULL) {
    info <- tibble(
        label =
            str_extract(colnames(as_tibble(normalised_gaps)), "\\w")
    ) %>%
        group_by(label) %>%
        summarise(count = dplyr::n()) %>%
        ungroup() %>%
        arrange(factor(label, levels = c("u", "t", "d"))) %>%
        mutate(count = cumsum(count))
    x_axis_labels <- c(1, pull(info, count))
    normalised_signal_tibble <- as.tibble(normalised_signal) %>%
        replace(is.na(.), 0) %>%
        rename_with(function(x) paste0("region_", colnames(.))) %>%
        pivot_longer(
            cols = starts_with("region_"), names_to = "position",
            values_to = "score", names_prefix = "region_"
        ) %>%
        mutate(
            coordinate = as.numeric(str_extract(position, "\\d+")),
            label = str_extract(position, "\\w"),
            coordinate =
                ifelse(label == "t",
                    coordinate + pull(filter(info, label == "u"), count),
                    ifelse(label == "d", coordinate + pull(filter(info, label == "t"), count),
                        coordinate
                    )
                )
        ) %>%
        group_by(coordinate) %>%
        summarise(signal_mean = mean(score))
    normalised_tibble <- as.tibble(normalised_gaps) %>%
        replace(is.na(.), 0) %>%
        rename_with(function(x) paste0("region_", colnames(.))) %>%
        pivot_longer(
            cols = starts_with("region_"), names_to = "position",
            values_to = "score", names_prefix = "region_"
        ) %>%
        mutate(
            coordinate = as.numeric(str_extract(position, "\\d+")),
            label = str_extract(position, "\\w"),
            coordinate =
                ifelse(label == "t",
                    coordinate + pull(filter(info, label == "u"), count),
                    ifelse(label == "d", coordinate + pull(filter(info, label == "t"), count),
                        coordinate
                    )
                )
        ) %>%
        group_by(coordinate) %>%
            summarise(mean = mean(score)) %>%
            ungroup() %>%
            inner_join(normalised_signal_tibble) %>%
            mutate(mean = mean / signal_mean) %>%
            select(-signal_mean) %>%
        mutate(r_mean = slide_dbl(mean, ~ mean(.x, na.rm = TRUE),
            .before = 2, .after = 2
        ), sample = label) %>%
        ungroup()
    annotation_plot <- ggplot(normalised_tibble, aes(
        x = coordinate, y = r_mean
    )) +
        geom_line() +
        annotate("rect",
            xmin = x_axis_labels[2], xmax = x_axis_labels[3], ymin = -Inf,
            ymax = Inf, alpha = 0.1, fill = "black"
        ) +
        scale_x_continuous(
            labels = c(extend[1], "TSS", "TES", extend[2]),
            breaks = x_axis_labels
        ) +
        labs(
            x = NULL, y = y_title
        ) +
        theme_cowplot() +
        grids(linetype = "dashed") +
        theme(
            strip.background = element_blank(),
            axis.line = element_blank(),
            panel.border = element_rect(colour = "black")
        )
    if (!is.null(path_for_output) > 0 && !is.null(prefix)) {
        ggsave(
            filename = paste0(path_for_output, prefix, ".pdf"), plot = annotation_plot,
            bg = "white", height = 10, width = 10
        )
    }
    return(annotation_plot)
}

# Plot distribution around TSS
plot_distribution_tss <- function(normalised_gaps, normalised_signal, extend = c(500, 1000),
                                  label = NULL,
                                  y_title = NULL,
                                  path_for_output = NULL,
                                  prefix = NULL) {
    info <- tibble(
        label =
            str_extract(colnames(as_tibble(normalised_gaps)), "\\w")
    ) %>%
        group_by(label) %>%
        summarise(count = dplyr::n()) %>%
        ungroup() %>%
        arrange(factor(label, levels = c("u", "d"))) %>%
        mutate(count = cumsum(count))
    x_axis_labels <- c(1, pull(info, count))
    normalised_signal_tibble <- as.tibble(normalised_signal) %>%
        replace(is.na(.), 0) %>%
        rename_with(function(x) paste0("region_", colnames(.))) %>%
        pivot_longer(
            cols = starts_with("region_"), names_to = "position",
            values_to = "score", names_prefix = "region_"
        ) %>%
        mutate(
            coordinate = as.numeric(str_extract(position, "\\d+")),
            label = str_extract(position, "\\w"),
            coordinate =
                ifelse(label == "d",
                    coordinate + pull(filter(info, label == "u"), count),
                    coordinate
                )
        ) %>%
        group_by(coordinate) %>%
        summarise(score_mean = mean(score))
    normalised_tibble <- as.tibble(normalised_gaps) %>%
        replace(is.na(.), 0) %>%
        rename_with(function(x) paste0("region_", colnames(.))) %>%
        pivot_longer(
            cols = starts_with("region_"), names_to = "position",
            values_to = "score", names_prefix = "region_"
        ) %>%
        mutate(
            coordinate = as.numeric(str_extract(position, "\\d+")),
            label = str_extract(position, "\\w"),
            coordinate =
                ifelse(label == "d",
                    coordinate + pull(filter(info, label == "u"), count),
                    coordinate
                )
        ) %>%
        group_by(coordinate) %>%
            summarise(mean = mean(score)) %>%
            inner_join(normalised_signal_tibble, by = c("coordinate" = "coordinate")) %>%
            ungroup() %>%
            mutate(mean = mean / score_mean) %>%
            select(-score_mean) %>%
            mutate(r_mean = slide_dbl(mean, ~ mean(.x, na.rm = TRUE),
            .before = 2, .after = 2
        ), sample = label)
    annotation_plot <- ggplot(normalised_tibble, aes(
        x = coordinate, y = r_mean
    )) +
        geom_line() +
        geom_vline(xintercept = x_axis_labels[2], colour = "black", lty = "dashed") +
        scale_x_continuous(
            labels = c(extend[1], "TSS", extend[2]),
            breaks = x_axis_labels
        ) +
        labs(
            x = NULL, y = y_title
        ) +
        theme_cowplot() +
        grids(linetype = "dashed") +
        theme(
            strip.background = element_blank(),
            axis.line = element_blank(),
            panel.border = element_rect(colour = "black")
        )
    if (!is.null(path_for_output) > 0 && !is.null(prefix)) {
        ggsave(
            filename = paste0(path_for_output, prefix, ".pdf"), plot = annotation_plot,
            bg = "white", height = 10, width = 10
        )
    }
    return(annotation_plot)
}


# Splitting gff into different features
split_gff <- function(path_to_gff, drop_seq_levels = NULL, convert_seq_levels = NULL) {
    transposable_elements <- import(path_to_gff) %>%
        dropSeqlevels(drop_seq_levels, "coarse") %>%
        filter(type == "transposable_element")
    txdb <- makeTxDbFromGFF(path_to_gff) %>%
        dropSeqlevels(drop_seq_levels, "coarse")
    if (!is.null(convert_seq_levels)) {
        seqlevels(transposable_elements) <- as.character(convert_seq_levels)
        seqlevels(txdb) <- as.character(convert_seq_levels)
    }
    genes <- genes(txdb)
    transposable_elements <- setdiff_ranges(transposable_elements, genes)
    exons <- exons(txdb)
    introns <- intronicParts(txdb)
    one_kb_upstream <- setdiff_ranges(promoters(genes,
        upstream = 1000,
        downstream = 0
    ), genes)
    one_kb_downstream <- genes
    start(one_kb_downstream[strand(one_kb_downstream) == "+"]) <-
        end(one_kb_downstream[strand(one_kb_downstream) == "+"]) + 1
    end(one_kb_downstream[strand(one_kb_downstream) == "+"]) <-
        end(one_kb_downstream[strand(one_kb_downstream) == "+"]) + 1000
    end(one_kb_downstream[strand(one_kb_downstream) == "-"]) <-
        start(one_kb_downstream[strand(one_kb_downstream) == "-"]) + 1
    start(one_kb_downstream[strand(one_kb_downstream) == "-"]) <-
        start(one_kb_downstream[strand(one_kb_downstream) == "-"]) - 1000
    one_kb_downstream <- setdiff_ranges(one_kb_downstream, genes)
    one_kb_downstream <- setdiff_ranges(one_kb_downstream, one_kb_upstream)
    one_kb_upstream <- setdiff_ranges(one_kb_upstream, one_kb_downstream)
    one_kb_downstream <- setdiff_ranges(one_kb_downstream, one_kb_upstream)
    output <- list(
        "genes" = genes,
        "exons" = exons,
        "introns" = introns,
        "transposable_elements" = transposable_elements,
        "one_kb_upstream" = one_kb_upstream,
        "one_kb_downstream" = one_kb_downstream
    )
    return(output)
}

# Annoating gaps to split features
annotate_gaps_to_features <- function(gaps_all, split_features, genome) {
    expected <- tibble(
        Feature = c(
            "Exons",
            "Introns",
            "Transposable Elements",
            "1kb Upstream",
            "1kb Downstream",
            "Intergenic"
        ),
        condition = "Expected",
        Length = c(
            sum(width(reduce_ranges(split_features$exons))),
            sum(width(reduce_ranges(split_features$introns))),
            sum(width(reduce_ranges(split_features$transposable_elements))),
            sum(width(reduce_ranges(split_features$one_kb_upstream))),
            sum(width(reduce_ranges(split_features$one_kb_downstream))), 0
        )
    ) %>%
        mutate(
            Length = ifelse(Length == 0, sum(width(genome)) - sum(Length), Length),
            total_length = sum(Length), Percentage = (Length / total_length) * 100
        )

    annotated_gaps_all <- c()
    ## Annotating all gaps (exons)
    annotated_gaps_all <- c(
        annotated_gaps_all,
        sum(width(reduce_ranges(intersect_ranges(
            gaps_all,
            split_features$exons
        ))))
    )
    gaps_all <- setdiff_ranges(gaps_all, split_features$exons)

    ## Annotating all gaps (introns)
    annotated_gaps_all <- c(
        annotated_gaps_all,
        sum(width(reduce_ranges(intersect_ranges(
            gaps_all,
            split_features$introns
        ))))
    )
    gaps_all <- setdiff_ranges(gaps_all, split_features$introns)

    ## Annotating all gaps (transposable_elements)
    annotated_gaps_all <- c(
        annotated_gaps_all,
        sum(width(reduce_ranges(intersect_ranges(
            gaps_all,
            split_features$transposable_elements
        ))))
    )
    gaps_all <- setdiff_ranges(
        gaps_all,
        split_features$transposable_elements
    )

    ## Annotating all gaps (one_kb_upstream)
    annotated_gaps_all <- c(
        annotated_gaps_all,
        sum(width(reduce_ranges(intersect_ranges(
            gaps_all,
            split_features$one_kb_upstream
        ))))
    )
    gaps_all <- setdiff_ranges(gaps_all, split_features$one_kb_upstream)

    ## Annotating all gaps (one_kb_downstream)
    annotated_gaps_all <- c(
        annotated_gaps_all,
        sum(width(reduce_ranges(intersect_ranges(
            gaps_all,
            split_features$one_kb_downstream
        ))))
    )
    gaps_all <- setdiff_ranges(gaps_all, split_features$one_kb_downstream)

    ## Annotating all gaps (Intergenic)
    annotated_gaps_all <- c(
        annotated_gaps_all,
        sum(width(reduce_ranges(gaps_all)))
    )

    annotated_gaps_all_df <- tibble(
        Feature = c(
            "Exons", "Introns", "Transposable Elements", "1kb Upstream",
            "1kb Downstream", "Intergenic"
        ),
        condition = "TF",
        Length = c(
            annotated_gaps_all
        )
    )
    annotated_gaps_all_df <- annotated_gaps_all_df %>%
        dplyr::summarise(total_length = sum(Length)) %>%
        bind_cols(annotated_gaps_all_df) %>%
        mutate(Percentage = (Length / total_length) * 100) %>%
        bind_rows(expected)
    return(annotated_gaps_all_df)
}

# Plotting annotated gaps
plot_annotated_gaps <- function(annotated_gaps, path_for_output = NULL, prefix = NULL, colour_palette ) {
    if (length(unique(annotated_gaps$context)) > 1) {
        output <- annotated_gaps %>%
            ggplot(aes(x = Percentage, y = condition, fill = Feature)) +
            geom_bar(stat = "identity", colour = "black") +
            # scale_fill_brewer(palette = "GnBu") +
            scale_fill_manual(values = colour_palette) +
            labs(y = NULL) +
            facet_rep_grid(rows = vars(context), repeat.tick.labels = TRUE) +
            theme_cowplot() +
            theme(axis.line.y = element_blank(), strip.background = element_blank())
    } else {
        output <- annotated_gaps %>%
            ggplot(aes(x = Percentage, y = condition, fill = Feature)) +
            geom_bar(stat = "identity", colour = "black") +
            # scale_fill_brewer(palette = "GnBu") +
            scale_fill_manual(values = colour_palette) +
            labs(y = NULL) +
            theme_cowplot() +
            theme(axis.line.y = element_blank(), strip.background = element_blank())
    }
    
    if (!is.null(path_for_output) > 0 && !is.null(prefix)) {
        ggsave(
            filename = paste0(path_for_output, prefix, ".pdf"), plot = output,
            bg = "white", height = 10, width = 10
        )
    }
    return(output)
}

# Function to simplify numbers
si_number <- function(x, digits = 2) {
    compress <- function(x, n) { # nolint
        signif(x * 10^(-n), digits)
    }
    case_when( # nolint
        x >= 1e6 ~ paste0(compress(x, 6), "M"),
        x >= 1000 ~ paste0(compress(x, 3), "k"),
        x >= 1 ~ as.character(compress(x, 0)),
        x >= 0.001 ~ paste0(compress(x, -3), "m"),
        x >= 1e-6 ~ paste0(compress(x, -6), "u")
    )
}

# Compute correlation between genes and gaps
compute_correlation <- function(gaps, genes, bin_width = 1e6) {
    binned_genes <- genes %>%
        as_tibble() %>%
        group_by(seqnames, bins = cut_width(start,
            width = bin_width,
            boundary = 0
        )) %>%
        summarise(n_genes = n(), .groups = "drop") %>%
        mutate(,
            n_genes = 100 * jitter(n_genes, amount = .5) / length(genes)
        )
    binned_gaps <- gaps %>%
        as_tibble() %>%
        group_by(seqnames, bins = cut_width(start,
            width = bin_width,
            boundary = 0
        )) %>%
        summarise(n_gaps = n(), .groups = "drop") %>%
        mutate(,
            n_gaps = 100 * jitter(n_gaps, amount = .5) / length(gaps)
        )

    binned_joined <- full_join(binned_genes, binned_gaps,
        by = c("seqnames", "bins")
    ) %>%
        na.omit()
    return(binned_joined)
}

# Plot correlation between genes and gaps
plot_gene_gap_correlation <- function(gene_gap_correlation, bin_width = 1e6, path_for_output = NULL, prefix = NULL, x_axis_title) {
    output <- ggplot(gene_gap_correlation, aes(x = n_genes, y = n_gaps)) +
        geom_smooth(method = "lm", se = FALSE, formula = y ~ x, size = 0.5) +
        geom_point(shape = 16, size = 1, alpha = 0.4) +
        scale_x_continuous(
            name = str_c(x_axis_title, si_number(bin_width), "bp"),
            expand = c(0, 0), limits = c(0, NA)
        ) +
        scale_y_continuous(
            name = str_c("% of TFs per ", si_number(bin_width), "bp"),
            expand = c(0, 0), limits = c(0, NA)
        ) +
        stat_regline_equation(
            label.y = 0.9 * max(gene_gap_correlation$n_gaps),
            label.x = 0.1 * max(gene_gap_correlation$n_genes)
        ) +
        stat_cor(
            label.y = 0.8 * max(gene_gap_correlation$n_gaps),
            label.x = 0.1 * max(gene_gap_correlation$n_genes)
        ) +
        facet_rep_grid(cols = vars(context), repeat.tick.labels = TRUE, scales = "free") +
        theme_cowplot() +
        theme(
            strip.background = element_blank(),
            axis.line = element_blank(),
            panel.border = element_rect(colour = "black", fill = NA, size = 1)
        )
    if (!is.null(path_for_output) > 0 && !is.null(prefix)) {
        ggsave(
            filename = paste0(path_for_output, prefix, ".pdf"), plot = output,
            bg = "white", height = 10, width = 10
        )
    }
}

# Load and convert narrowpeak to gff3
read_narrow_peak <- function(path_to_file, convert_to_gff3 = FALSE,
                             path_for_output = NULL, prefix = NULL) {
    narrow_peak <- read_delim(path_to_file,
        col_names = c(
            "seqnames", "start", "end", "name", "score", "strand",
            "signalValue", "pValue", "qValue", "peak"
        )
    ) %>%
        makeGRangesFromDataFrame(keep.extra.columns = TRUE)
    if (convert_to_gff3) {
        narrow_peak <- narrow_peak %>%
            select(name, score)
    }
    if (convert_to_gff3 & !is.null(path_for_output) & !is.null(prefix)) {
        export(
            narrow_peak,
            paste0(
                path_for_output, prefix,
                ".gff3"
            )
        )
    }
    return(narrow_peak)
}

# Plotting overlap as barchart
overlap_data <- function(data1, data2, data1_title, path_for_output = NULL, prefix = NULL) {
    overlap_tibble <- tibble()
    for (i in seq_len(length(data1))) {
        area1 <- setdiff_ranges(data1[[i]], data2) %>%
            width() %>%
            sum()
        area2 <- setdiff_ranges(data2, data1[[i]]) %>%
            width() %>%
            sum()
        cross_area <- intersect_ranges(data1[[i]], data2) %>%
            width() %>%
            sum()
        overlap_tibble <- bind_rows(overlap_tibble, tibble(
            label = names(split_gaps)[i],
            category = c(data1_title, "overlap"),
            length = c(area1, cross_area), total = sum(length),
            percentage = (length / total) * 100
        ))
    }
    output <- list()
    if (length(data1) == 1) {
        overlap_tibble <- overlap_tibble %>%
            mutate(label = "Merged")
    }
    output[[1]] <- overlap_tibble
    output[[2]] <- ggplot(overlap_tibble, aes(x = percentage, y = label, fill = category)) +
        geom_bar(stat = "identity", colour = "black") +
        geom_text(aes(label = signif(percentage, 3)),
            color = "black",
            size = 5, position = position_stack(vjust = 0.5)
        ) +
        labs(x = "Percentage", y = "Context") +
        theme_cowplot() +
        grids(linetype = "dashed") +
        theme(panel.border = element_rect(colour = "black"))
    if (!is.null(path_for_output) & !is.null(prefix)) {
        ggsave(filename = paste0(path_for_output, prefix, ".pdf"), plot = output[[2]], bg = "white", height = 10, width = 10)
    }
    return(output)
}
