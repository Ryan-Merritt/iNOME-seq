rm(list=ls())
source("~/iNOME/Scripts/R/libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
#### Methylation correlation for leaf ####
## Loading Fasta file
fasta <- import("~/work/ecotypes/annotations/Col0/TAIR10_chr_all.fasta")[1:5]
names(fasta) <- paste0("Chr",seq(1,5,1))
chromosome_lengths <- width(fasta)
## Computing methylation correlation
files <- dir(path = "/Volumes/PhD_one/NOME_Seq/data/CX_reports/leaf/GpC/")
contexts <- c("CG","CHG","CHH")
leaf_correlation_GC <- list()
leaf_correlation_non_GC <- list()
for(i in files){
  temp_CX <- read_delim(paste0("/Volumes/PhD_one/NOME_Seq/data/CX_reports/leaf/GpC/",i),delim = "\t",
                        col_names = c("seqnames","start","end","width","strand",
                                      "readsM","readsN","context","trinucleotide_context","sequence")) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    keepSeqlevels(paste0("Chr",c(seq(1,5,1))),"coarse")
  temp_GC <- temp_CX[temp_CX$sequence == "GC"]
  temp_CX <- temp_CX[temp_CX$sequence != "GC"]
  leaf_correlation_GC[[i]] <- list()
  leaf_correlation_non_GC[[i]] <- list()
  for(j in 1:3){
    leaf_correlation_GC[[i]][[j]] <- temp_GC %>%
      computeMethylationDataSpatialCorrelation(methylationData = .,context = contexts[j],distances = c(50,100,200,500,1000,2000,5000)) %>%
      tibble(position = seq(1,7,1),correlation = .,context = paste0("G",contexts[j]),sample = gsub("_CX_report.txt","",i))
    leaf_correlation_non_GC[[i]][[j]] <- temp_CX %>%
      computeMethylationDataSpatialCorrelation(methylationData = .,context = contexts[j],distances = c(50,100,200,500,1000,2000,5000)) %>%
      tibble(position = seq(1,7,1),correlation = .,context = paste0("G",contexts[j]),sample = gsub("_CX_report.txt","",i))
  }
  leaf_correlation_GC[[i]] <- bind_rows(leaf_correlation_GC[[i]])
  leaf_correlation_non_GC[[i]] <- bind_rows(leaf_correlation_non_GC[[i]])
}
leaf_correlation_GC <- bind_rows(leaf_correlation_GC)
leaf_correlation_non_GC <- bind_rows(leaf_correlation_non_GC)
# Plotting methylation correlation
leaf_correlation_GC <- leaf_correlation_GC %>%
  mutate(time = gsub(".*_","",sample),sample = gsub("_.*","",sample))
leaf_correlation_GC$time <- factor(leaf_correlation_GC$time,levels = c("6h","12h","24h"))

leaf_correlation_non_GC <- leaf_correlation_non_GC %>%
  mutate(time = gsub(".*_","",sample),sample = gsub("_.*","",sample))
leaf_correlation_non_GC$time <- factor(leaf_correlation_non_GC$time,levels = c("6h","12h","24h"))

colours <- c('#8b0000', '#a56252', '#a9a9a9', '#8e65b1', '#00008b')
ggplot(leaf_correlation_GC,aes(x = position,y = correlation,colour = sample)) +
  geom_line() +
  facet_grid(cols = vars(context),rows = vars(time)) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(labels = c(50,100,200,500,1000,2000,5000),breaks = seq(1,7,1)) +
  scale_y_continuous(labels = seq(0,1,0.1),breaks = seq(0,1,0.1)) +
  xlab("Distances between Cytosine bases (log(bp))") +
  geom_hline(yintercept = -Inf,size = 1) +
  theme_cowplot() +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,vjust = 0.5),axis.text = element_text(size = 8))
ggsave("~/iNOME/Figures/Supplementary/leaf_correlation_GC_Col0.pdf")

ggplot(leaf_correlation_non_GC,aes(x = position,y = correlation,colour = sample)) +
  geom_line() +
  facet_grid(cols = vars(context),rows = vars(time)) +
  scale_colour_manual(values = colours) +
  scale_x_continuous(labels = c(50,100,200,500,1000,2000,5000),breaks = seq(1,7,1)) +
  scale_y_continuous(labels = seq(0,1,0.1),breaks = seq(0,1,0.1)) +
  xlab("Distances between Cytosine bases (log(bp))") +
  geom_hline(yintercept = -Inf,size = 1) +
  theme_cowplot() +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,vjust = 0.5),axis.text = element_text(size = 8))
ggsave("~/iNOME/Figures/Supplementary/leaf_correlation_non_GC_Col0.pdf")




