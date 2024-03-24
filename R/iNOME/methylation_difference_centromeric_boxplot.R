rm(list=ls())
source("~/iNOME/Scripts/R/libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
# Creating Centromeric Regions --------------------------------------------
centromeric_regions <- tibble(seqnames = c("Chr1","Chr2","Chr3","Chr4","Chr4","Chr5"),
                              start = c(12500000,2300000,12800000,1620000,2780000,10680000),
                              end = c(17050000,6300000,14800000,2280000,5804000,14000000)) %>%
  makeGRangesFromDataFrame()
# Loading Genome ----------------------------------------------------------
## Importing genome from fasta file
file = "~/work/ecotypes/annotations/Col0/TAIR10_chr_all.fasta"
genome = import(file,format = "fasta")[1:5]
chromosome_lengths = width(genome)
names(genome) <- sub(" .*", "", names(genome))
# Loading CX reports from plant lines -------------------------------------
## Loading CX reports
path <- "/Volumes/PhD_one/NOME_Seq/Julia/cx_lines/"
files <- dir(path = path)
MCViP_lines <- lapply(paste0(path,files), function(x){
  x <- readBismark(x)
  x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  x <- sort(x)
  x
})
names(MCViP_lines) <- gsub(".R1.T_bismark_bt2_pe.deduplicated.CX_report.txt","",files)
# Extracting GC and Non-GC from plant lines -------------------------------
MCViP_GC <- lapply(MCViP_lines,function(x){
  x <- extract_seq(x,fasta = genome,lengths = chromosome_lengths,result = "GC") %>%
    filter(sequence == "GC")
  x
})
MCViP_Non_GC <- lapply(MCViP_lines,function(x){
  x <- extract_seq(x,fasta = genome,lengths = chromosome_lengths,result = "non_GC") %>%
    filter(sequence != "GC")
  x
})
# Computing methylation profiles ------------------------------------------
## Computing methylation profiles GC
window <- 750000
contexts <- c("CG","CHG","CHH")
MCViP_GC_profile <- vector("list",length(MCViP_GC))
names(MCViP_GC_profile) <- names(MCViP_GC)
for(i in 1:length(MCViP_GC_profile)){
  regions <- getWholeChromosomes(MCViP_GC[[i]])
  MCViP_GC_profile[[i]] <- vector("list",length(regions))
  names(MCViP_GC_profile[[i]]) <- seqlevels(regions)
  for(j in 1:length(MCViP_GC_profile[[i]])){
    MCViP_GC_profile[[i]][[j]] <- vector("list",length(contexts))
    names(MCViP_GC_profile[[i]][[j]]) <- contexts
    for(k in 1:length(MCViP_GC_profile[[i]][[j]])){
      MCViP_GC_profile[[i]][[j]][[k]] <- 
        computeMethylationProfile(MCViP_GC[[i]],region = regions[j],
                                  windowSize = window,context = contexts[k])
    }
    MCViP_GC_profile[[i]][[j]] <- bind_ranges(MCViP_GC_profile[[i]][[j]])
  }
  MCViP_GC_profile[[i]] <- bind_ranges(MCViP_GC_profile[[i]])
}
## Computing methylation profiles non-GC
MCViP_Non_GC_profile <- vector("list",length(MCViP_Non_GC))
names(MCViP_Non_GC_profile) <- names(MCViP_Non_GC)
for(i in 1:length(MCViP_Non_GC_profile)){
  regions <- getWholeChromosomes(MCViP_Non_GC[[i]])
  MCViP_Non_GC_profile[[i]] <- vector("list",length(regions))
  names(MCViP_Non_GC_profile[[i]]) <- seqlevels(regions)
  for(j in 1:length(MCViP_Non_GC_profile[[i]])){
    MCViP_Non_GC_profile[[i]][[j]] <- vector("list",length(contexts))
    names(MCViP_Non_GC_profile[[i]][[j]]) <- contexts
    for(k in 1:length(MCViP_Non_GC_profile[[i]][[j]])){
      MCViP_Non_GC_profile[[i]][[j]][[k]] <- 
        computeMethylationProfile(MCViP_Non_GC[[i]],region = regions[j],
                                  windowSize = window,context = contexts[k])
    }
    MCViP_Non_GC_profile[[i]][[j]] <- bind_ranges(MCViP_Non_GC_profile[[i]][[j]])
  }
  MCViP_Non_GC_profile[[i]] <- bind_ranges(MCViP_Non_GC_profile[[i]])
}
# Extracting Centromeric and non-centromeric cytosines --------------------
centromeric_cytosines_GC <- list()
centromeric_cytosines_non_GC <- list()
for(i in 1:length(MCViP_GC_profile)){
  centromeric_cytosines_GC[[i]] <- subsetByOverlaps(MCViP_GC_profile[[i]],centromeric_regions) %>%
    as_tibble() %>%
    select(seqnames,context,Proportion) %>%
    mutate(sequence = "GC",line = names(MCViP_GC_profile)[i])
  centromeric_cytosines_non_GC[[i]] <- subsetByOverlaps(MCViP_Non_GC_profile[[i]],centromeric_regions) %>%
    as_tibble() %>%
    select(seqnames,context,Proportion) %>%
    mutate(sequence = "Non-GC",line = names(MCViP_Non_GC_profile)[i])
}
centromeric_cytosines_GC <- bind_rows(centromeric_cytosines_GC) %>%
  mutate(location = "Centromeric")
centromeric_cytosines_non_GC <- bind_rows(centromeric_cytosines_non_GC) %>%
  mutate(location = "Centromeric")

non_centromeric_cytosines_GC <- list()
non_centromeric_cytosines_non_GC <- list()
for(i in 1:length(MCViP_GC_profile)){
  non_centromeric_cytosines_GC[[i]] <- subsetByOverlaps(MCViP_GC_profile[[i]],centromeric_regions,invert = T) %>%
    as_tibble() %>%
    select(seqnames,context,Proportion) %>%
    mutate(sequence = "GC",line = names(MCViP_GC_profile)[i])
  non_centromeric_cytosines_non_GC[[i]] <- subsetByOverlaps(MCViP_Non_GC_profile[[i]],centromeric_regions,invert = T) %>%
    as_tibble() %>%
    select(seqnames,context,Proportion) %>%
    mutate(sequence = "Non-GC",line = names(MCViP_Non_GC_profile)[i])
}
non_centromeric_cytosines_GC <- bind_rows(non_centromeric_cytosines_GC) %>%
  mutate(location = "Non-Centromeric")
non_centromeric_cytosines_non_GC <- bind_rows(non_centromeric_cytosines_non_GC) %>%
  mutate(location = "Non-Centromeric")

all_cytosines <- bind_rows(centromeric_cytosines_GC,centromeric_cytosines_non_GC,non_centromeric_cytosines_GC,non_centromeric_cytosines_non_GC)

# Plotting boxplots for each comparison
comparisons <- tibble(reference = "LCol_0E",treatment = c("L10aE","L1aE","L2aE","L5aE"))
colours <- c('#8b0000', '#a56252', '#a9a9a9', '#8e65b1')
for(i in 1:nrow(comparisons)){
  reference_line <- comparisons %>%
    slice(i) %>%
    pull(reference)
  treatment_line <- comparisons %>%
    slice(i) %>%
    pull(treatment)
  x <- all_cytosines %>%
    filter(line == reference_line | line == treatment_line,seqnames == "Chr1") %>%
    rename(Line = line,Location = location)
  x$Line <- factor(x$Line,levels = c(reference_line,treatment_line))
  stats <- x %>%
    group_by(Location,context,sequence) %>%
    wilcox_test(Proportion ~ Line) %>%
    add_xy_position(x = "Location") %>%
    #mutate(y.position = rep(1.1,12)) %>%
    add_significance()
  y <- ggboxplot(x, x = "Location", y = "Proportion", fill = "Line",
            scales = "free",outlier.shape = NA,size = 0.1) +
    facet_rep_grid(cols = vars(context),rows = vars(sequence),scales = "free",repeat.tick.labels = T) +
    stat_pvalue_manual(stats,label = "p.signif",tip.length = 0,bracket.size = 0.1) +
    scale_fill_manual(values = c('#00008b',colours[[i]])) +
    theme_cowplot() +
    theme(strip.background = element_blank(),strip.text.y = element_text(angle = 0),axis.title.x = element_blank(),
          axis.line = element_line(size = 0.1),axis.ticks = element_line(size = 0.1))
  ggsave(plot = y,filename = paste0("~/iNOME/Figures/Main/Figure1/Parts/",reference_line,".",treatment_line,".pdf"),width = 12)
}

# Plotting all comparisons
x <- all_cytosines %>%
  filter(line != "LCol_0") %>%
  rename(Line = line,Location = location)
x$Line <- factor(x$Line,levels = c("L10aE","L1aE","L2aE","L5aE","LCol_0E"))
colours <- c('#8b0000', '#a56252', '#a9a9a9', '#8e65b1', '#00008b')

stats <- x %>%
  group_by(Location,context,sequence) %>%
  wilcox_test(Proportion ~ Line) %>%
  add_xy_position(x = "Location")

ggboxplot(x, x = "Location", y = "Proportion", fill = "Line",
          scales = "free",outlier.shape = NA,size = 0.1) +
  facet_rep_grid(cols = vars(context),rows = vars(sequence),scales = "free",repeat.tick.labels = T) +
  stat_pvalue_manual(stats,label = "p.adj.signif",tip.length = 0,bracket.size = 0.1) +
  scale_fill_manual(values = colours) +
  theme_cowplot() +
  theme(strip.background = element_blank(),strip.text.y = element_text(angle = 0),axis.title.x = element_blank(),
        axis.line = element_line(size = 0.1),axis.ticks = element_line(size = 0.1))
ggsave("~/iNOME/Figures/Main/Figure1/Parts/all_boxplot.pdf",width = 12)
