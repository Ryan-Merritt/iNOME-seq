rm(list=ls())
source("~/iNOME/Scripts/R/libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
#### Methylation profile for root ####
## Loading Fasta file
fasta <- import("~/work/ecotypes/annotations/Col0/TAIR10_chr_all.fasta")[1:5]
names(fasta) <- paste0("Chr",seq(1,5,1))
chromosome_lengths <- width(fasta)
## Computing methylation profile
files <- dir(path = "/Volumes/PhD_one/NOME_Seq/data/CX_reports/root/GpC/")
contexts <- c("CG","CHG","CHH")
root_profile_GC <- list()
root_profile_non_GC <- list()
for(i in files){
  temp_CX <- read_delim(paste0("/Volumes/PhD_one/NOME_Seq/data/CX_reports/root/GpC/",i),delim = "\t",
                        col_names = c("seqnames","start","end","width","strand",
                                      "context","readsM","readsN","trinucleotide_context","sequence")) %>%
    mutate(end = start) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    keepSeqlevels(paste0("Chr",c(seq(1,5,1))),"coarse")
  regions <- getWholeChromosomes(temp_CX)
  root_profile_GC[[i]] <- list()
  root_profile_non_GC[[i]] <- list()
  for(j in 1:3){
    root_profile_GC[[i]][[j]] <- list()
    root_profile_non_GC[[i]][[j]] <- list()
    for(k in 1:5){
      root_profile_non_GC[[i]][[j]][[k]] <- temp_CX %>%
        filter(context == contexts[j],sequence != "GC") %>%
        computeMethylationProfile(methylationData = .,context = contexts[j],region = regions[k],window = 100000) %>%
        as_tibble() %>%
        summarise(chromosome = seqnames,position = (start + end)/2, proportion = Proportion,context = contexts[j],sample = gsub("_CX_report.txt","",i))
      root_profile_GC[[i]][[j]][[k]] <- temp_CX %>%
        filter(context == contexts[j],sequence == "GC") %>%
        computeMethylationProfile(methylationData = .,context = contexts[j],region = regions[k],window = 100000) %>%
        as_tibble() %>%
        summarise(chromosome = seqnames,position = (start + end)/2, proportion = Proportion,context = paste0("G",contexts[j]),sample = gsub("_CX_report.txt","",i))
    }
    root_profile_GC[[i]][[j]] <- bind_rows(root_profile_GC[[i]][[j]])
    root_profile_non_GC[[i]][[j]] <- bind_rows(root_profile_non_GC[[i]][[j]])
  }
  root_profile_GC[[i]] <- bind_rows(root_profile_GC[[i]])
  root_profile_non_GC[[i]] <- bind_rows(root_profile_non_GC[[i]])
}
root_profile_GC <- bind_rows(root_profile_GC)
root_profile_non_GC <- bind_rows(root_profile_non_GC)
# Plotting methylation profile
colours <- c('#8b0000', '#a56252', '#a9a9a9', '#8e65b1', '#00008b')
root_profile_GC <- root_profile_GC %>%
  mutate(time = gsub(".*_","",sample),sample = gsub("_.*","",sample))
root_profile_GC$time <- factor(root_profile_GC$time,levels = c("6h","12h","24h"))

root_profile_non_GC <- root_profile_non_GC %>%
  mutate(time = gsub(".*_","",sample),sample = gsub("_.*","",sample))
root_profile_non_GC$time <- factor(root_profile_non_GC$time,levels = c("6h","12h","24h"))

times <- c("6h","12h","24h")
for(i in times){
  root_profile_GC %>%
    filter(time == i) %>%
    ggplot(aes(x = position/1e6,y = proportion,colour = sample)) +
    geom_line() +
    facet_grid(cols = vars(chromosome),rows = vars(context),scales = "free_x") +
    scale_colour_manual(values = colours) +
    scale_x_continuous(labels = c(seq(0,30,5)),breaks = c(seq(0,30,5))) +
    ylim(0,1) + 
    ylab("Methylation Proportion") +
    xlab("Genomic Coordinate (Mb)") +
    geom_hline(yintercept = -Inf,size = 1) +
    theme_cowplot() +
    theme(strip.background = element_blank(),axis.text = element_text(size = 8))
  ggsave(paste0("~/iNOME/Figures/Supplementary/root_profile_",i,"_GC_Col0.pdf"))
}

for(i in times){
  root_profile_non_GC %>%
    filter(time == i) %>%
    ggplot(aes(x = position/1e6,y = proportion,colour = sample)) +
    geom_line() +
    facet_grid(cols = vars(chromosome),rows = vars(context),scales = "free_x") +
    scale_colour_manual(values = colours) +
    scale_x_continuous(labels = c(seq(0,30,5)),breaks = c(seq(0,30,5))) +
    ylim(0,1) + 
    ylab("Methylation Proportion") +
    xlab("Genomic Coordinate (Mb)") +
    geom_hline(yintercept = -Inf,size = 1) +
    theme_cowplot() +
    theme(strip.background = element_blank(),axis.text = element_text(size = 8))
  ggsave(paste0("~/iNOME/Figures/Supplementary/root_profile_",i,"_non_GC_Col0.pdf"))
}




