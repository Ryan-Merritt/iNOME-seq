rm(list=ls())
source("~/iNOME/Scripts/R/libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
#### Methylation coverage for root ####
## Loading Fasta file
fasta <- import("~/work/ecotypes/annotations/Col0/TAIR10_chr_all.fasta")[1:5]
names(fasta) <- paste0("Chr",seq(1,5,1))
chromosome_lengths <- width(fasta)
## Computing methylation coverage
files <- dir(path = "/Volumes/PhD_one/NOME_Seq/data/CX_reports/root/GpC/")
contexts <- c("CG","CHG","CHH")
root_coverage_GC <- list()
root_coverage_non_GC <- list()
for(i in files){
  temp_CX <- read_delim(paste0("/Volumes/PhD_one/NOME_Seq/data/CX_reports/root/GpC/",i),delim = "\t",
                        col_names = c("seqnames","start","end","width","strand","context",
                                      "readsM","readsN","trinucleotide_context","sequence")) %>%
    mutate(end = start) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    keepSeqlevels(paste0("Chr",c(seq(1,5,1))),"coarse") 
  root_coverage_GC[[i]] <- list()
  root_coverage_non_GC[[i]] <- list()
  for(j in 1:3){
    root_coverage_non_GC[[i]][[j]] <- temp_CX %>%
      filter(context == contexts[j],sequence != "GC") %>%
      computeMethylationDataCoverage(methylationData = .,context = contexts[j],breaks = seq(1,15,1)) %>%
      tibble(position = seq(1,15,1),coverage = .,context = contexts[j],sample = gsub("_CX_report.txt","",i))
    root_coverage_GC[[i]][[j]] <- temp_CX %>%
      filter(context == contexts[j],sequence == "GC") %>%
      computeMethylationDataCoverage(methylationData = .,context = contexts[j],breaks = seq(1,15,1)) %>%
      tibble(position = seq(1,15,1),coverage = .,context = paste0("G",contexts[j]),sample = gsub("_CX_report.txt","",i))
  }
  root_coverage_GC[[i]] <- bind_rows(root_coverage_GC[[i]])
  root_coverage_non_GC[[i]] <- bind_rows(root_coverage_non_GC[[i]])
}
root_coverage_non_GC <- bind_rows(root_coverage_non_GC)
root_coverage_GC <- bind_rows(root_coverage_GC)
# Plotting methylation coverage
colours <- c('#8b0000', '#a56252', '#a9a9a9', '#8e65b1', '#00008b')
root_coverage_GC <- root_coverage_GC %>%
  mutate(time = gsub(".*_","",sample),sample = gsub("_.*","",sample))
root_coverage_GC$time <- factor(root_coverage_GC$time,levels = c("6h","12h","24h"))

root_coverage_non_GC <- root_coverage_non_GC %>%
  mutate(time = gsub(".*_","",sample),sample = gsub("_.*","",sample))
root_coverage_non_GC$time <- factor(root_coverage_non_GC$time,levels = c("6h","12h","24h"))

ggplot(root_coverage_GC,aes(x = position,y = coverage,colour = sample)) +
  geom_line() +
  facet_grid(cols = vars(context),rows = vars(time)) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(labels = c(seq(1,15,1)),breaks = c(seq(1,15,1))) +
  xlab("Minimum number of reads") +
  geom_hline(yintercept = -Inf,size = 1) +
  theme_cowplot() +
  theme(strip.background = element_blank(),axis.text = element_text(size = 8))
ggsave("~/iNOME/Figures/Supplementary/root_coverage_GC_Col0.pdf")

colours <- c('#8b0000', '#a56252', '#a9a9a9', '#8e65b1', '#00008b')
root_coverage_non_GC %>%
  ggplot(aes(x = position,y = coverage,colour = sample)) +
  geom_line() +
  facet_grid(cols = vars(context),rows = vars(time)) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(labels = c(seq(1,15,1)),breaks = c(seq(1,15,1))) +
  xlab("Minimum number of reads") +
  geom_hline(yintercept = -Inf,size = 1) +
  theme_cowplot() +
  theme(strip.background = element_blank(),axis.text = element_text(size = 8))
ggsave("~/iNOME/Figures/Supplementary/root_coverage_non_GC_Col0.pdf")





