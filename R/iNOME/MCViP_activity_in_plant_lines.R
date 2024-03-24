source("~/iNOME/Scripts/R/libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
#### Checking MCViP activity in plant lines
# Loading Genome ----------------------------------------------------------
## Importing genome from fasta file\
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
  x <- extract_seq(x,fasta = genome,lengths = chromosome_lengths) %>%
    filter(sequence == "GC")
  x
})
MCViP_Non_GC <- lapply(MCViP_lines,function(x){
  x <- extract_seq(x,fasta = genome,lengths = chromosome_lengths) %>%
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
  }
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
  }
}
# Plotting MCViP activity in plant lines (3) ----------------------------------
## Creating empty data frame
MCViP_df <- data.frame(Chromosome=factor(),
                       Proportion = numeric(),
                       Coordinate = numeric(),
                       Context = factor(),
                       Line = factor(),
                       Type = factor(),
                       stringsAsFactors=FALSE)
## Combining data frames to allow for easy plotting
for(i in 1:length(MCViP_GC_profile)){
  for(j in 1:length(MCViP_GC_profile[[i]])){
    for(k in 1:length(MCViP_GC_profile[[i]][[j]])){
      MCViP_df <- rbind(MCViP_df,data.frame(Chromosome = seqnames(MCViP_GC_profile[[i]][[j]][[k]]),
                                            Coordinate = (start(MCViP_GC_profile[[i]][[j]][[k]]) + end(MCViP_GC_profile[[i]][[j]][[k]]))/2,
                                            Proportion = MCViP_GC_profile[[i]][[j]][[k]]$Proportion,
                                            Context = MCViP_GC_profile[[i]][[j]][[k]]$context,
                                            Line = names(MCViP_GC_profile)[i],
                                            Type = "GC"))
      MCViP_df <- rbind(MCViP_df,data.frame(Chromosome = seqnames(MCViP_Non_GC_profile[[i]][[j]][[k]]),
                                            Coordinate = (start(MCViP_Non_GC_profile[[i]][[j]][[k]]) + end(MCViP_Non_GC_profile[[i]][[j]][[k]]))/2,
                                            Proportion = MCViP_Non_GC_profile[[i]][[j]][[k]]$Proportion,
                                            Context = MCViP_Non_GC_profile[[i]][[j]][[k]]$context,
                                            Line = names(MCViP_Non_GC_profile)[i],
                                            Type = "non-GC"))
    }
  }
}
## Removing Irrelevent data
colours <- c('#8b0000', '#a56252', '#a9a9a9', '#8e65b1', '#00008b')
MCViP_df %>%
  filter(Chromosome == "Chr1", Line != "LCol_0") %>%
  ggplot(aes(x = Coordinate/1000000,y = Proportion,colour = Line)) +
  geom_line() +
  facet_rep_grid(cols = vars(Context),rows = vars(Type),repeat.tick.labels = T) +
  scale_color_manual(values = c(colours)) +
  xlab("Genomic Coordinate (Mb)") +
  ylab("Methylation Proportion") +
  theme_cowplot() +
  theme(strip.background = element_blank())
ggsave("~/iNOME/Figures/main/Figure1/Parts/MCViP_activity_500_smooth.pdf")
