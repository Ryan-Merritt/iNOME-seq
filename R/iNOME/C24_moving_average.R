source("~/iNOME/Scripts/R/libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
fasta <- import("~/work/Julia/ecotype_comparisons_hub/C24/Annotation/C24_ragtag.fasta")[1:5]
chromosome_lengths <- width(fasta)
# Computing Control methylation -------------------------------------------
files <- dir("~/iNOME/Data/cx_lines/C24/CX/",pattern = "0h",full.names = T)
control_methylation <- read_delim(files,delim = "\t",
                                  col_names = c("seqnames","start","strand","readsM","readsN","context","trinucleotide_context")) %>%
  mutate(end = start,score = readsM/(readsM + readsN)) %>%
  filter(is.finite(score)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
split_control_methylation <- split(control_methylation,control_methylation$context)
for(i in 1:3){
  split_control_methylation[[i]] <- split_control_methylation[[i]] %>%
    select(score)
  seqlengths(split_control_methylation[[i]]) <- end(getWholeChromosomes(control_methylation))
  export(split_control_methylation[[i]],paste0("~/iNOME/Data/cx_lines/C24/control_methylation/C24_control_",names(split_control_methylation)[i],".bw"))
}

# Computing moving proportion (all) ---------------------------------------
# extracting GpC
files <- dir("~/iNOME/Data/cx_lines/C24/CX/")
for(i in files){
  x <- read_delim(paste0("~/iNOME/Data/cx_lines/C24/CX/",i),delim = "\t",
                  col_names = c("seqnames","start","strand","readsM","readsN","context","trinucleotide_context")) %>%
    mutate(end = start) %>%
    filter(seqnames != "ChrM",seqnames != "ChrC") %>%
    extract_seq(.,fasta,chromosome_lengths)
  write.table(x,paste0("~/iNOME/Data/cx_lines/C24/GpC/",i),append = F,quote = F,
              sep = "\t",row.names = F,col.names = F)
}
# Loading GpC data
files <- dir("~/iNOME/Data/cx_lines/C24/GpC/")
GpC <- list()
for(i in files){
  GpC[[i]] <- read_delim(paste0("~/iNOME/Data/cx_lines/C24/GpC/",i),delim = "\t",
                  col_names = c("seqnames","start","end","width","strand","readsM","readsN","context","trinucleotide_context","sequence")) %>%
    filter(sequence == "GC") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  GpC[[i]] <- split(GpC[[i]],GpC[[i]]$context)
} 
names(GpC) <- c("0h","12h","6h")
# Calculating moving proportion
pb <- txtProgressBar(0,length(GpC),style = 3)
GC_prop <- vector("list",length(GpC))
names(GC_prop) <- names(GpC)
for(i in 1:length(GpC)){
  GC_prop[[i]] <- movingproportion(GpC[[i]],500)
  GC_prop[[i]] <- GC_prop[[i]][seq(1,nrow(GC_prop[[i]]),10),]
  colnames(GC_prop[[i]]) <- c("chromosome","start","end","score")
  GC_prop[[i]] <- GRanges(seqnames = GC_prop[[i]]$chromosome,
                              ranges = IRanges(GC_prop[[i]]$start,
                                               (GC_prop[[i]]$end + 9)),
                              score = GC_prop[[i]]$score)
  seqlengths(GC_prop[[i]]) <- end(getWholeChromosomes(GC_prop[[i]]))
  setTxtProgressBar(pb,i)
}
save(GC_prop,file = paste0("~/iNOME/Data/RData/C24_500_ma_all.rdata"))

control <- GC_prop[[1]]
treatment <- GC_prop[[2]]

treatment$score <- treatment$score - control$score
treatment <- treatment[!is.na(treatment$score)]
treatment <- treatment[is.finite(treatment$score)]

export(treatment,paste0("~/iNOME/Data/cx_lines/C24/bw/C24_500_ma_all.bw"),"bw")
# Calculate moving proportion (split) -------------------------------------
# Loading GpC data
files <- dir("~/iNOME/Data/cx_lines/C24/GpC/")
GpC <- list()
for(i in files){
  GpC[[i]] <- read_delim(paste0("~/iNOME/Data/cx_lines/C24/GpC/",i),delim = "\t",
                         col_names = c("seqnames","start","end","width","strand","readsM","readsN","context","trinucleotide_context","sequence")) %>%
    filter(sequence == "GC") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
} 
names(GpC) <- c("0h","12h","6h")
# Calculating moving proportion
pb <- txtProgressBar(0,length(GpC),style = 3)
GC_prop <- vector("list",length(GpC))
names(GC_prop) <- names(GpC)
for(i in 1:length(GpC)){
  x <- GpC[[i]] %>%
    filter(context == "CG")
  GC_prop[[i]] <- movingproportion(x,500)
  GC_prop[[i]] <- GC_prop[[i]][seq(1,nrow(GC_prop[[i]]),10),]
  colnames(GC_prop[[i]]) <- c("chromosome","start","end","score")
  GC_prop[[i]] <- GRanges(seqnames = GC_prop[[i]]$chromosome,
                          ranges = IRanges(GC_prop[[i]]$start,
                                           (GC_prop[[i]]$end + 9)),
                          score = GC_prop[[i]]$score)
  seqlengths(GC_prop[[i]]) <- end(getWholeChromosomes(GC_prop[[i]]))
  setTxtProgressBar(pb,i)
}
save(GC_prop,file = paste0("~/iNOME/Data/RData/C24_500_ma_all.rdata"))
