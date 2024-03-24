rm(list=ls())
#####################
#### Preparation ####
#####################
####-####
# Installing required packages --------------------------------------------
if(!requireNamespace("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
bioc.packages <- c("rtracklayer","EnrichedHeatmap","GenomicFeatures","DMRcaller",
                   "BSgenome","compEpiTools","DESeq2","EnhancedVolcano",
                   "plyranges","tximport","BSgenome.Athaliana.TAIR.TAIR9")
bioconductor.check <- lapply(bioc.packages,FUN = function(x){
  if(!require(x,character.only = TRUE)){
    BiocManager::install(x,character.only = TRUE)
    library(x,character.only = T)
  }
})

packages <- c("devtools", "RColorBrewer", "ggplot2", "ggpubr", "cowplot", 
              "circlize", "scales", "VennDiagram", "pheatmap", "matrixStats", 
              "arrangements", "ggforce", "ggalluvial","tidyverse","ggsci")
packages.check <- lapply(packages,FUN = function(x){
  if(!require(x,character.only = TRUE)){
    install.packages(x,character.only = TRUE,quiet = T)
    library(x,character.only = T)
  }
})
rm(bioc.packages,bioconductor.check,packages.check,packages)
# Functions ---------------------------------------------------------------
# Colour function for heatmaps
brewer.ramp = function(start, end, pal) {
  palette = brewer.pal(n = brewer.pal.info[pal, "maxcolors"],  name = pal)
  sequence = seq(start, end, (end - start)/(brewer.pal.info[pal, "maxcolors"] - 1))
  ramp = colorRamp2(sequence, palette)
  return(ramp)
}
col_fun_me3 = brewer.ramp(0, 3, "Reds")

# Moving average
.stopIfNotAll <- function(exprs, errorMsg) {
  for(expr in exprs){
    if (! expr)
      stop(errorMsg, call. = FALSE)
  }
}

movingAverage <- function(minPos, maxPos, pos, val, weights=1, windowSizeHalf=150, normalize = FALSE, kernelFunction="triangular", lambda=0.5) {
  
  .stopIfNotAll(c(length(pos)==length(val)), "pos and val vectors need to have the same length")
  
  if(length(weights) < length(val)){
    weights <- rep(weights, length.out=length(val))
  } else if(length(weights) > length(val)){
    weights <- weights[1:length(val)]
  }
  
  # Filter out NAs
  keepIndexes <- which(!is.na(pos) & !is.na(val) & !is.na(weights))
  pos <- pos[keepIndexes]
  weights <- weights[keepIndexes]  
  val <- val[keepIndexes]
  
  #set the values
  rawVector <- rep(0, maxPos - minPos + 2*windowSizeHalf + 1)
  rawVector[(pos - minPos + windowSizeHalf + 1)] <-weights*val
  
  normVector <- rep(0, maxPos - minPos + 2*windowSizeHalf + 1)
  normVector[(pos - minPos + windowSizeHalf + 1)] <-weights
  
  
  
  # Define the (triangular) kernel.
  if(kernelFunction =="uniform"){
    kernel <- rep(1, times=(2*windowSizeHalf + 1))  
  } else if(kernelFunction =="triangular"){
    kernel <- 1 - abs(-windowSizeHalf:windowSizeHalf)/windowSizeHalf 
  } else if(kernelFunction =="gaussian"){
    kernel <- .gaussianKernel(lambda, - windowSizeHalf:windowSizeHalf)
  } else if(kernelFunction =="epanechnicov"){
    kernel <- .epanechnicovKernel(-windowSizeHalf:windowSizeHalf, windowSizeHalf)
  } else{
    stop(paste("Unknown kernel function: ", kernelFunction, ". 
               It should be one of \"uniform\", \"triangular\", \"gaussian\", \"epanechnicov\"",sep=""))
  }
  
  kernel <- kernel / sum(kernel)
  
  if(windowSizeHalf >= 1){
    smoothedVector <- RcppRoll::roll_sum(rawVector, length(kernel), weights = kernel, normalize = normalize) / RcppRoll::roll_sum(normVector, length(kernel), weights = kernel, normalize = normalize)
  } else{
    smoothedVector <- rawVector 
  }
  
  return(smoothedVector);
}

# Moving Proportion
movingproportion <- function(file,windowsize){
  chromosomes <- seqlevels(file)
  proportion <- vector("list",length(chromosomes))
  for(i in 1:length(chromosomes)){
    chrData <- file[seqnames(file)==chromosomes[[i]]]
    region <- GRanges(seqnames = unique(seqnames(chrData)),
                      ranges = IRanges(start = min(start(chrData)),
                                       max(end(chrData))))
    windowSizeHalf <- floor((windowsize -1)/2)
    methylR1 <- movingAverage(start(region),end(region),start(chrData),
                              chrData$readsM,windowSizeHalf = windowSizeHalf)
    totalReads1 <- movingAverage(start(region),end(region),start(chrData),
                                 chrData$readsN,windowSizeHalf = windowSizeHalf)
    movingProp <- methylR1 / totalReads1
    if(length(is.na(movingProp))>0){
      movingProp[is.na(movingProp)] <- 0
    }
    positions <- start(region) : end(region)
    proportion[[i]] <- data.frame(chromosomes[[i]],positions,positions,
                                  movingProp)
  }
  result <- do.call("rbind",proportion)
  return(result)
}
# Extract genes from txdb
genes_txdb <- GenomicFeatures::genes
# Loading CX_reports
readBismark2 <- function(file){
  .stopIfNotAll(c(!is.null(file), is.character(file), file.exists(file)), 
                " file does not exist")
  cat("Reading file: ", file, "\n", sep = "")
  cx <- scan(file = file, what = list(character(), integer(), 
                                      character(), integer(), integer(), character(), character()), 
             sep = "\t", skip = 0)
  chrs <- unique(cx[[1]])
  chrs <- chrs[order(chrs)]
  dat <- data.frame(chr = factor(cx[[1]], levels = chrs), pos = as.integer(cx[[2]]), 
                    strand = factor(cx[[3]], levels = c("+", "-")), context = factor(cx[[6]], 
                                                                                     levels = c("CG", "CHG", "CHH")), trinucleotide_context = cx[[7]], 
                    readsM = as.integer(cx[[4]]), readsN = as.integer(cx[[5]]))
  dat <- GRanges(seqnames = dat$chr, ranges = IRanges(start = dat$pos, 
                                                      end = dat$pos), strand = dat$strand, context = dat$context, 
                 readsM = dat$readsM, readsN = dat$readsN, trinucleotide_context = dat$trinucleotide_context)
  cat("Finished reading file: ", file, "\n", sep = "")
  return(dat)
}

geom_cytoband <- function(cent = centromeres, hchrom = heterochromatin, chr_col = chr) {
  cent   <- rename(cent,   {{ chr_col }} := chr)
  hchrom <- rename(hchrom, {{ chr_col }} := chr)
  list(geom_vline(data = cent, aes(xintercept = position/1e+6), linetype="dashed", color = "black"),
       geom_rect(data = hchrom, aes(x = NULL, y = NULL, xmin = start/1e+6, xmax = end/1e+6, ymin = -Inf, ymax = Inf), 
                 color = NA, fill = "black", alpha = 0.1))
}

facet_chr <- function(chr = chr, rows = NULL, scales = "free", space = "free_x", breaks = 10) {
  list(facet_grid(rows = vars({{ rows }}), cols = vars({{ chr }}), scales = scales, space = space),
       scale_x_continuous(expand = c(0, 0), limits = c(0, NA), breaks = seq(0, 1000, breaks)), 
       xlab("Position across the genome (Mbp)"))
}
.remove_bad_cytosines <- function(x,lengths){
  for(i in 1:5){
    x[[i]] <- x[[i]] %>%
      filter(!end > lengths[i],start > 0)
  }
  return(x)
}

.Extract_sequence <- function(x,fasta){
  temp.storage <- getSeq(fasta,x)
  x <- x %>%
    mutate(sequence = as.character(temp.storage))
  return(x)
}

convert_cx_data <- function(x,set.context = NULL){
  seqlevels(x) <- paste0("Chr",seqlevels(x))
  if(!is.null(set.context)){
    x$context <- set.context
  }
  return(x)
}

tss_focus <- function(x,down,up){
  x <- as_tibble(x) %>%
    mutate(new.start = ifelse(strand == "+",start - up,end + up),
           new.end = ifelse(strand == "+",start + down,end - down),
           start = ifelse(strand == "+",new.start,new.end),
           end = ifelse(strand == "+",new.end,new.start)) %>%
    select(-new.start,-new.end) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  return(x)
}

combine_GRanges <- function(x,y,z){
  x <- as_tibble(x)
  y <- as_tibble(y)
  o <- inner_join(x,y,by = z)
  return(o)
}

unique_overlap_ranges <- function(gff,regions){
  average_DMR_proportion <- find_overlaps(gff,regions) %>%
    group_by(ID) %>%
    summarise(score = mean(proportion))
  gff <- inner_join(as_tibble(gff),as_tibble(average_DMR_proportion),by = "ID") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  return(gff)
}

extract_GC <- function(methylation.data,fasta,lengths){
  output <- as_tibble(methylation.data) %>%
    mutate(start = ifelse(strand == "+",start - 1,start),end = ifelse(strand == "-",end + 1,end)) %>%
    group_split(seqnames) %>%
    .remove_bad_cytosines(.,lengths) %>%
    bind_rows() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    .Extract_sequence(.,fasta) %>%
    filter(sequence == "GC") %>%
    mutate(start = ifelse(strand == "+",start + 1,start),end = ifelse(strand == "-",end - 1,end))
  return(output)
}

####-####
##########################
#### Loading data (1,2) ####
##########################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs <- root_12h_CHH_ME_CE$`0.3`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE)
# Loading ATAC peaks -------------------------------------------------------
## Loading ATAC-seq peaks
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak")
ATAC_peaks <- keepSeqlevels(ATAC_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
ATAC_peaks <- GRanges(seqnames(ATAC_peaks),IRanges(start(ATAC_peaks),end(ATAC_peaks)))
ATAC_peaks <- unique(ATAC_peaks)
# Loading DNase peaks ------------------------------------------------------
## Loading DNase-seq peaks
DNase_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/DNase_q005_peaks.narrowPeak")
seqlevels(DNase_peaks) <- paste0("Chr",seqlevels(DNase_peaks))
DNase_peaks <- keepSeqlevels(DNase_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
DNase_peaks <- GRanges(seqnames(DNase_peaks),IRanges(start(DNase_peaks),end(DNase_peaks)))
DNase_peaks <- unique(DNase_peaks)
# Loading GFF -------------------------------------------------------------
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Splitting DMRs ----------------------------------------------------------
## Extracting common DMRs
cDMRs <- intersect(DMRs,ATAC_peaks)
cDMRs <- intersect(cDMRs,DNase_peaks)
## Extracting unique DMRs
sDMRs <- setdiff(DMRs,ATAC_peaks)
sDMRs <- setdiff(sDMRs,DNase_peaks)
## Saving backup
cDMRs_backup <- cDMRs
sDMRs_backup <- sDMRs
# Loading bigwigs ---------------------------------------------------------
## Loading ATAC data
ATAC <- import("/Volumes/PhD_one/NOME_Seq/Col_features/bw/ATAC_col_10.bw")
ATAC <- keepSeqlevels(ATAC,c("Chr1","Chr2","Chr3","Chr4","Chr5"),pruning.mode = "coarse")
ATAC$type <- "ATAC"
## Loading DNase data
DNase <- import("/Volumes/PhD_one/NOME_Seq/Col_features/bw/DNase_col_10.bw")
DNase <- keepSeqlevels(DNase,c("Chr1","Chr2","Chr3","Chr4","Chr5"),pruning.mode = "coarse")
DNase$type <- "DNase"
## Loading H3K27me3
H3K27me3 <- import("/Volumes/PhD_one/NOME_Seq/Col_features/bw/H3K27me3_col_10.bw")
H3K27me3 <- keepSeqlevels(H3K27me3,c("Chr1","Chr2","Chr3","Chr4","Chr5"),pruning.mode = "coarse")
H3K27me3$type <- "H3K27me3"
## Loading H3K27ac
H3K27ac <- import("/Volumes/PhD_one/NOME_Seq/Col_features/bw/h3k27ac_col_10.bw")
H3K27ac <- keepSeqlevels(H3K27ac,c("Chr1","Chr2","Chr3","Chr4","Chr5"),pruning.mode = "coarse")
H3K27ac$type <- "H3K27ac"
## Combining data
feature_scores <- GRangesList(ATAC,DNase,H3K27me3,H3K27ac)
## Cleaning global environment
rm(ATAC,DNase,H3K27me3,H3K27ac)
####-####
##########################################
#### Normalising data around DMRs (2) ####
##########################################
####-####
# Normalising data around cDMRs -------------------------------------------
## Changing cDMRs to midpoint
cDMRs_backup$midpoint <- (start(cDMRs_backup) + end(cDMRs_backup))/2
start(cDMRs_backup) <- cDMRs_backup$midpoint
end(cDMRs_backup) <- cDMRs_backup$midpoint
## Normalising data around cDMRs
pb <- txtProgressBar(0,length(feature_scores),style = 3)
normalised_cDMRs <- vector("list",length(feature_scores))
for(i in 1:length(feature_scores)){
  normalised_cDMRs[[i]] <- normalizeToMatrix(feature_scores[[i]],cDMRs_backup,
                                             value_column = "score",
                                             mean_mode = "absolute",
                                             extend = 2000,w = 10,
                                             background = NA,smooth = F)
  setTxtProgressBar(pb,i)
}

# Normalising data around sDMRs -------------------------------------------
## Changing sDMRs to midpoint
sDMRs_backup$midpoint <- (start(sDMRs_backup) + end(sDMRs_backup))/2
start(sDMRs_backup) <- sDMRs_backup$midpoint
end(sDMRs_backup) <- sDMRs_backup$midpoint
## Normalising data around sDMRs
pb <- txtProgressBar(0,length(feature_scores),style = 3)
normalised_sDMRs <- vector("list",length(feature_scores))
for(i in 1:length(feature_scores)){
  normalised_sDMRs[[i]] <- normalizeToMatrix(feature_scores[[i]],sDMRs_backup,
                                             value_column = "score",
                                             mean_mode = "absolute",
                                             extend = 2000,w = 10,
                                             background = NA,smooth = F)
  setTxtProgressBar(pb,i)
}
####-####
####################################################
#### Checking MCViP activity in plant lines (3) ####
####################################################
####-####
# Loading Genome ----------------------------------------------------------
## Importing genome from fasta file\
file = "~/work/ecotypes/annotations/Col0/TAIR10_chr_all.fasta"
genome = import(file,format = "fasta")[1:5]
names(genome) <- sub(" .*", "", names(genome))
## Cleaning global environment
rm(file)
# Loading CX reports from plant lines -------------------------------------
## Loading CX reports
path <- "~/work/Julia/cx_lines/"
files <- dir(path = path)
MCViP_lines <- lapply(paste0(path,files), function(x){
  x <- readBismark(x)
  x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  x <- sort(x)
  x
})
names(MCViP_lines) <- gsub(".R1.T_bismark_bt2_pe.deduplicated.CX_report.txt","",files)
## Cleaning global environment
rm(files,path)
# Extracting GC and Non-GC from plant lines -------------------------------
MCViP_GC <- lapply(MCViP_lines,extractGC,genome = genome,contexts = "ALL")
MCViP_GC <- lapply(MCViP_GC,unlist)
MCViP_Non_GC <- lapply(MCViP_lines,extractNonGC,genome = genome,contexts = "ALL")
MCViP_Non_GC <- lapply(MCViP_Non_GC,unlist)
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
####-####
##################################################
#### Methylation profiles of root vs leaf (4) ####
##################################################
####-####
# Loading Genome ----------------------------------------------------------
## Importing genome from fasta file\
file = "~/work/ecotypes/annotations/Col0/TAIR10_chr_all.fasta"
genome = import(file,format = "fasta")[1:5]
names(genome) <- sub(" .*", "", names(genome))
## Cleaning global environment
rm(file)
# Loading data ------------------------------------------------------------
## Loading root CX report
if(!file.exists("~/work/RData/root_GC_unsplit.RData")){
  path <- "~/work/Julia/cx_root_vs_leaf/root/"
  files <- dir(path = path)
  root_meth <- lapply(paste0(path,files),function(x){
    x <- readBismark(x)
    seqlevels(x) <- paste0("Chr",seqlevels(x))
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
    x <- sort(x)
  })
  names(root_meth) <- gsub(pattern = "_R.R1.RUN0076.L3._CX_report.txt",replacement = "",x = files)
}

## Loading Leaf CX report
if(!file.exists("~/work/RData/leaf_GC_unsplit.RData")){
  path <- "~/work/Julia/cx_root_vs_leaf/leaf/"
  files <- dir(path = path)
  leaf_meth <- lapply(paste0(path,files),function(x){
    x <- readBismark(x)
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
    x <- sort(x)
  })
  names(leaf_meth) <- gsub(pattern = "_CX_report.txt",replacement = "",x = files)
}
# Extracting GC -----------------------------------------------------------
## Extracting root GC
if(!file.exists("~/work/RData/root_GC_unsplit.RData")){
  root_meth <- lapply(root_meth,extractGC,genome = genome,contexts = "ALL")
  root_GC_unsplit <- lapply(root_meth,unlist)
  save(root_GC_unsplit,file = "~/work/RData/root_GC_unsplit.RData")
} else{
  load("~/work/RData/root_GC_unsplit.RData")
}
## Extracting leaf GC
if(!file.exists("~/work/RData/leaf_GC_unsplit.RData")){
  leaf_meth <- lapply(leaf_meth,extractGC,genome = genome,contexts = "ALL")
  leaf_GC_unsplit <- lapply(leaf_meth,unlist)
  save(leaf_GC_unsplit,file = "~/work/RData/leaf_GC_unsplit.RData")
} else{
  load("~/work/RData/leaf_GC_unsplit.RData")
}
# Computing methylation profiles ------------------------------------------
## Computing root methylation profiles
window <- 500000
contexts <- c("CG","CHG","CHH")
root_profile <- vector("list",length(root_GC_unsplit))
names(root_profile) <- names(root_GC_unsplit)
for(i in 1:length(root_profile)){
  regions <- getWholeChromosomes(root_GC_unsplit[[i]])
  root_profile[[i]] <- vector("list",length(regions))
  names(root_profile[[i]]) <- seqlevels(regions)
  for(j in 1:length(root_profile[[i]])){
    root_profile[[i]][[j]] <- vector("list",length(contexts))
    names(root_profile[[i]][[j]]) <- contexts
    for(k in 1:length(root_profile[[i]][[j]])){
      root_profile[[i]][[j]][[k]] <- 
        computeMethylationProfile(root_GC_unsplit[[i]],region = regions[j],
                                  windowSize = window,context = contexts[k])
    }
  }
}
## Computing leaf methylation profiles
leaf_profile <- vector("list",length(leaf_GC_unsplit))
names(leaf_profile) <- names(leaf_GC_unsplit)
for(i in 1:length(leaf_profile)){
  regions <- getWholeChromosomes(leaf_GC_unsplit[[i]])
  leaf_profile[[i]] <- vector("list",length(regions))
  names(leaf_profile[[i]]) <- seqlevels(regions)
  for(j in 1:length(leaf_profile[[i]])){
    leaf_profile[[i]][[j]] <- vector("list",length(contexts))
    names(leaf_profile[[i]][[j]]) <- contexts
    for(k in 1:length(leaf_profile[[i]][[j]])){
      leaf_profile[[i]][[j]][[k]] <- 
        computeMethylationProfile(leaf_GC_unsplit[[i]],region = regions[j],
                                  windowSize = window,context = contexts[k])
    }
  }
}
## Cleaning global environment
rm(contexts,i,j,k,window,regions)
####-####
########################################
#### Overlap with other methods (5) ####
########################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CpG_ME_CE_12h.RData")
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/leaf/CpG_ME_CE_12h.RData")

load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHG_ME_CE_12h.RData")
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/leaf/CHG_ME_CE_12h.RData")

load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/leaf/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs_root_CpG <- root_12h_CpG_ME_CE$`0.3`$`700`
DMRs_leaf_CpG <- leaf_12h_CpG_ME_CE$`0.15`$`700`

DMRs_root_CHG <- root_12h_CHG_ME_CE$`0.3`$`700`
DMRs_leaf_CHG <- leaf_12h_CHG_ME_CE$`0.15`$`700`

DMRs_root_CHH <- root_12h_CHH_ME_CE$`0.3`$`700`
DMRs_leaf_CHH <- leaf_12h_CHH_ME_CE$`0.15`$`700`
## Cleaning global environment
rm("leaf_12h_CHG_ME_CE", "leaf_12h_CHH_ME_CE", 
   "leaf_12h_CpG_ME_CE", "root_12h_CHG_ME_CE", "root_12h_CHH_ME_CE", 
   "root_12h_CpG_ME_CE")
# Loading data from other methods -----------------------------------------
## Loading ATAC-seq data (Root)
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak")
ATAC_peaks <- keepSeqlevels(ATAC_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
ATAC_peaks <- GRanges(seqnames(ATAC_peaks),IRanges(start(ATAC_peaks),end(ATAC_peaks)))
ATAC_peaks <- unique(ATAC_peaks)
# Overlapping with other methods data -------------------------------------
## Storing full width of accessible regions detected by methods (Root)
total_width_root_CpG <- data.frame(Method = c("iNOME","ATAC"),
                                   Length = c(sum(width(DMRs_root_CpG)),sum(width(ATAC_peaks))))
total_width_root_CpG$label <- round(total_width_root_CpG$Length/1000000,1)

total_width_root_CHG <- data.frame(Method = c("iNOME"," ATAC"),
                                   Length = c(sum(width(DMRs_root_CHG)),sum(width(ATAC_peaks))))
total_width_root_CHG$label <- round(total_width_root_CHG$Length/1000000,1)

total_width_root_CHH <- data.frame(Method = c("iNOME","ATAC"),
                                   Length = c(sum(width(DMRs_root_CHH)),sum(width(ATAC_peaks))))
total_width_root_CHH$label <- round(total_width_root_CHH$Length/1000000,1)

## Overlapping DMRs with other methods (Root)
DMRs_ATAC_root_CpG <- sum(width(intersect(DMRs_root_CpG,ATAC_peaks)))

DMRs_ATAC_root_CHG <- sum(width(intersect(DMRs_root_CHG,ATAC_peaks)))

DMRs_ATAC_root_CHH <- sum(width(intersect(DMRs_root_CHH,ATAC_peaks)))
####-####
#################################################
#### Annotating root and leaf DMRs (6) DNase ####
#################################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/leaf/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs_root <- root_12h_CHH_ME_CE$`0.3`$`700`
DMRs_leaf <- leaf_12h_CHH_ME_CE$`0.15`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE,leaf_12h_CHH_ME_CE)
# Loading GFF -------------------------------------------------------------
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Loading GFF (Transposable elements) -------------------------------------
TE <- import("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
TE <- keepSeqlevels(TE,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
TE <- TE[TE$type=="transposable_element"]

# Preparation -------------------------------------------------------------
## Creating data frame to store the ammount of overlaps
root_overlap <- data.frame(Type = character(),Length = integer())
leaf_overlap <- data.frame(Type = character(),Length = integer())
## Storing total length of accessible regions detected
root_overlap <- rbind(root_overlap,data.frame(Type = "Total",Length = sum(width(DMRs_root))))
leaf_overlap <- rbind(leaf_overlap,data.frame(Type = "Total",Length = sum(width(DMRs_leaf))))
# Genes -------------------------------------------------------------------
## Storing length of DMRs that overlap with genes
root_overlap <- rbind(root_overlap,data.frame(Type = "Genes",Length = sum(width(intersect(DMRs_root,genes(GFF),ignore.strand = T)))))
leaf_overlap <- rbind(leaf_overlap,data.frame(Type = "Genes",Length = sum(width(intersect(DMRs_leaf,genes(GFF),ignore.strand = T)))))
## Removing regions that overlapped with genes from total DMRs
DMRs_root <- setdiff(DMRs_root,genes(GFF),ignore.strand = T)
DMRs_leaf <- setdiff(DMRs_leaf,genes(GFF),ignore.strand = T)
# Transposable elements ---------------------------------------------------
## Storing length of DMRs that overlap with Transposable elements
root_overlap <- rbind(root_overlap,data.frame(Type = "Transposable Elements",Length = sum(width(intersect(DMRs_root,TE,ignore.strand = T)))))
leaf_overlap <- rbind(leaf_overlap,data.frame(Type = "Transposable Elements",Length = sum(width(intersect(DMRs_leaf,TE,ignore.strand = T)))))
## Removing regions that overlapped with Transposable elements from total DMRs
DMRs_root <- setdiff(DMRs_root,TE,ignore.strand = T)
DMRs_leaf <- setdiff(DMRs_leaf,TE,ignore.strand = T)
# 1kb upstream ------------------------------------------------------------
## Storing length of DMRs that are less that 1kb upstream of genes
root_overlap <- rbind(root_overlap,data.frame(Type = "1kb upstream",Length = sum(width(intersect(DMRs_root,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)))))
leaf_overlap <- rbind(leaf_overlap,data.frame(Type = "1kb upstream",Length = sum(width(intersect(DMRs_leaf,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)))))
## Removing regions that are less than 1kb away of genes from total DMRs
DMRs_root <- setdiff(DMRs_root,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)
DMRs_leaf <- setdiff(DMRs_leaf,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)
# 1kb downstream ----------------------------------------------------------
## Extracting downstream from GFF
one_kb_downstream <- genes(GFF)
start(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"])
end(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"]) + 1000
end(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"])
start(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"]) - 1000
## Storing length of DMRs that are less that 1kb downstream of genes
root_overlap <- rbind(root_overlap,data.frame(Type = "1kb downstream",Length = sum(width(intersect(DMRs_root,one_kb_downstream,ignore.strand = T)))))
leaf_overlap <- rbind(leaf_overlap,data.frame(Type = "1kb downstream",Length = sum(width(intersect(DMRs_leaf,one_kb_downstream,ignore.strand = T)))))
## Removing regions that are less than 1kb away of genes from total DMRs
DMRs_root <- setdiff(DMRs_root,one_kb_downstream,ignore.strand = T)
DMRs_leaf <- setdiff(DMRs_leaf,one_kb_downstream,ignore.strand = T)
# Intergenic --------------------------------------------------------------
## Storing length of DMRs that are less that 1kb downstream of genes
root_overlap <- rbind(root_overlap,data.frame(Type = "Intergenic",Length = sum(width(DMRs_root))))
leaf_overlap <- rbind(leaf_overlap,data.frame(Type = "Intergenic",Length = sum(width(DMRs_leaf))))

####-####
##########################################
#### Comparing DMRs to expression (7) ####
##########################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs_root <- root_12h_CHH_ME_CE$`0.3`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE)
# Loading RNA-seq expression ----------------------------------------------
## Loading RNA data and removing irrelevent data
RNA_root <- import("/Volumes/PhD_one/NOME_Seq/Col_features/bw/RNA_root_edited.bw")
RNA_root <- keepSeqlevels(RNA_root,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Calculating DMR score and average expression ----------------------------
## Calculating DMR score
DMRs_root$score <- DMRs_root$proportion2 - DMRs_root$proportion1
## Removing irrelevent data
DMRs_root <- GRanges(seqnames(DMRs_root),IRanges(start(DMRs_root),end(DMRs_root)),score = DMRs_root$score)
## Computing average expression for each DMR
DMRs_root$expression <- 0
pb <- txtProgressBar(0,length(DMRs_root),style = 3)
for(i in 1:length(DMRs_root)){
  DMRs_root[i]$expression <- mean(subsetByOverlaps(RNA_root,DMRs_root[i])$score)
  setTxtProgressBar(pb,i)
}
df <- data.frame(DMRs = DMRs_root$score,RNA = DMRs_root$expression)
ggplot(df,aes(x = RNA,y = DMRs)) + 
  geom_point() +
  xlim(0,2) +
  ylim(0,1)
####-####
#############################################
#### iNOME_seq around genes metaplot (8) ####
#############################################
####-####
# Loading data ------------------------------------------------------------
## Loading CX reports
if(!file.exists("~/work/RData/root_GC_unsplit.RData")){
  path <- "~/work/Julia/cx_root_vs_leaf/root/"
  files <- dir(path = path)
  root_meth <- lapply(paste0(path,files),function(x){
    x <- readBismark(x)
    seqlevels(x) <- paste0("Chr",seqlevels(x))
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
    x <- sort(x)
  })
  names(root_meth) <- gsub(pattern = "_R.R1.RUN0076.L3._CX_report.txt",replacement = "",x = files)
}
## Loading GC data
if(!file.exists("~/work/RData/root_GC_split.RData")){
  root_meth <- lapply(root_meth,extractGC,genome = genome,contexts = "ALL")
  save(root_GC_unsplit,file = "~/work/RData/root_GC_unsplit.RData")
} else{
  load("~/work/RData/root_GC_split.RData")
}
# Loading GFF -------------------------------------------------------------
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Calculating moving proportion -------------------------------------------
if(!file.exists("~/work/RData/root_GC_prop.RData")){
  GC_prop <- vector("list",length(root_GC))
  names(GC_prop) <- names(root_GC)
  for(i in 1:length(root_GC)){
    message(names(root_GC)[i])
    pb <- txtProgressBar(0,length(root_GC[[i]]),style = 3)
    GC_prop[[i]] <-vector("list",length(root_GC[[i]]))
    names(GC_prop[[i]]) <- c("CG","CHG","CHH")
    for(j in 1:length(root_GC[[i]])){
      GC_prop[[i]][[j]] <- movingproportion(root_GC[[i]][[j]],100)
      GC_prop[[i]][[j]] <- GC_prop[[i]][[j]][seq(1,nrow(GC_prop[[i]][[j]]),10),]
      colnames(GC_prop[[i]][[j]]) <- c("chromosome","start","end","score")
      GC_prop[[i]][[j]] <- GRanges(seqnames = GC_prop[[i]][[j]]$chromosome,
                                   ranges = IRanges(GC_prop[[i]][[j]]$start,
                                                    (GC_prop[[i]][[j]]$end + 9)),
                                   score = GC_prop[[i]][[j]]$score,
                                   minus_score = (1-GC_prop[[i]][[j]]$score))
      seqlengths(GC_prop[[i]][[j]]) <- end(getWholeChromosomes(GC_prop[[i]][[j]]))
      setTxtProgressBar(pb,j)
    }
    message("\n")
  }
  save(GC_prop,file = "~/work/RData/root_GC_prop.RData")
} else{
  ## Loading RData
  load("~/work/RData/root_GC_prop.RData")
}
# Normalizing data around features ----------------------------------------
GC_control <- GC_prop[grep(pattern = "CE",x = names(GC_prop))]
GC_treatment <- GC_prop[grep(pattern = "ME",x = names(GC_prop))]
for(i in 1:3){
  for(j in 1:3){
    GC_treatment[[i]][[j]]$score <- GC_treatment[[i]][[j]]$score - GC_control[[i]][[j]]$score
  }
}
GC_mat <- vector("list",length(GC_treatment))
names(GC_mat) <- names(GC_treatment)
for(i in 1:length(GC_treatment)){
  message(names(GC_treatment)[i])
  GC_mat[[i]] <- vector("list",length(GC_treatment[[i]]))
  names(GC_mat[[i]]) <- names(GC_treatment[[i]])
  pb <- txtProgressBar(0,length(GC_treatment[[i]]),style = 3)
  for(j in 1:length(GC_treatment[[i]])){
    GC_mat[[i]][[j]] <- normalizeToMatrix(GC_treatment[[i]][[j]],genes(GFF),
                                          value_column = "score",mean_mode = "absolute",
                                          extend = 1000, w = 10, smooth = F)
    GC_mat[[i]][[j]] <- as.data.frame(colMeans(GC_mat[[i]][[j]],na.rm = T))
    GC_mat[[i]][[j]] <- cbind(GC_mat[[i]][[j]],seq(10,3330,10),names(GC_treatment[i]),names(GC_treatment[[i]])[j])
    colnames(GC_mat[[i]][[j]]) <- c("Score","Position","Sample","Context")
    setTxtProgressBar(pb,j)
  }
}
####-####
#################################################
#### Annotating common and specific DMRs (9) ####
#################################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs <- root_12h_CHH_ME_CE$`0.3`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE)
# Loading ATAC peaks -------------------------------------------------------
## Loading ATAC-seq peaks
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak")
ATAC_peaks <- keepSeqlevels(ATAC_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
ATAC_peaks <- GRanges(seqnames(ATAC_peaks),IRanges(start(ATAC_peaks),end(ATAC_peaks)))
ATAC_peaks <- unique(ATAC_peaks)
# Loading DNase peaks ------------------------------------------------------
## Loading DNase-seq peaks
DNase_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/DNase_q005_peaks.narrowPeak")
seqlevels(DNase_peaks) <- paste0("Chr",seqlevels(DNase_peaks))
DNase_peaks <- keepSeqlevels(DNase_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
DNase_peaks <- GRanges(seqnames(DNase_peaks),IRanges(start(DNase_peaks),end(DNase_peaks)))
DNase_peaks <- unique(DNase_peaks)
# Loading GFF -------------------------------------------------------------
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Loading GFF (Transposable elements) -------------------------------------
TE <- import("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
TE <- keepSeqlevels(TE,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
TE <- TE[TE$type=="transposable_element"]
# Splitting DMRs ----------------------------------------------------------
## Extracting common DMRs
cDMRs <- intersect(DMRs,ATAC_peaks)
cDMRs <- intersect(cDMRs,DNase_peaks)
## Extracting unique DMRs
sDMRs <- setdiff(DMRs,ATAC_peaks)
sDMRs <- setdiff(sDMRs,DNase_peaks)
## Saving backup
cDMRs_backup <- cDMRs
sDMRs_backup <- sDMRs


# Preparation -------------------------------------------------------------
## Creating data frame to store the ammount of overlaps
cDMRs_overlap <- data.frame(Type = character(),Length = integer())
sDMRs_overlap <- data.frame(Type = character(),Length = integer())
## Storing total length of accessible regions detected
cDMRs_overlap <- rbind(cDMRs_overlap,data.frame(Type = "Total",Length = sum(width(cDMRs))))
sDMRs_overlap <- rbind(sDMRs_overlap,data.frame(Type = "Total",Length = sum(width(sDMRs))))
# Genes -------------------------------------------------------------------
## Storing length of DMRs that overlap with genes
cDMRs_overlap <- rbind(cDMRs_overlap,data.frame(Type = "Genes",Length = sum(width(intersect(cDMRs,genes(GFF),ignore.strand = T)))))
sDMRs_overlap <- rbind(sDMRs_overlap,data.frame(Type = "Genes",Length = sum(width(intersect(sDMRs,genes(GFF),ignore.strand = T)))))
## Removing regions that overlapped with genes from total DMRs
cDMRs <- setdiff(cDMRs,genes(GFF),ignore.strand = T)
sDMRs <- setdiff(sDMRs,genes(GFF),ignore.strand = T)
# Transposable elements ---------------------------------------------------
## Storing length of DMRs that overlap with Transposable elements
cDMRs_overlap <- rbind(cDMRs_overlap,data.frame(Type = "Transposable Elements",Length = sum(width(intersect(cDMRs,TE,ignore.strand = T)))))
sDMRs_overlap <- rbind(sDMRs_overlap,data.frame(Type = "Transposable Elements",Length = sum(width(intersect(sDMRs,TE,ignore.strand = T)))))
## Removing regions that overlapped with Transposable elements from total DMRs
cDMRs <- setdiff(cDMRs,TE,ignore.strand = T)
sDMRs <- setdiff(sDMRs,TE,ignore.strand = T)
# 1kb upstream ------------------------------------------------------------
## Storing length of DMRs that are less that 1kb upstream of genes
cDMRs_overlap <- rbind(cDMRs_overlap,data.frame(Type = "1kb upstream",Length = sum(width(intersect(cDMRs,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)))))
sDMRs_overlap <- rbind(sDMRs_overlap,data.frame(Type = "1kb upstream",Length = sum(width(intersect(sDMRs,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)))))
## Removing regions that are less than 1kb away of genes from total DMRs
cDMRs <- setdiff(cDMRs,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)
sDMRs <- setdiff(sDMRs,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)
# 1kb downstream ----------------------------------------------------------
## Extracting downstream from GFF
one_kb_downstream <- genes(GFF)
start(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"])
end(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"]) + 1000
end(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"])
start(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"]) - 1000
## Storing length of DMRs that are less that 1kb downstream of genes
cDMRs_overlap <- rbind(cDMRs_overlap,data.frame(Type = "1kb downstream",Length = sum(width(intersect(cDMRs,one_kb_downstream,ignore.strand = T)))))
sDMRs_overlap <- rbind(sDMRs_overlap,data.frame(Type = "1kb downstream",Length = sum(width(intersect(sDMRs,one_kb_downstream,ignore.strand = T)))))
## Removing regions that are less than 1kb away of genes from total DMRs
cDMRs <- setdiff(cDMRs,one_kb_downstream,ignore.strand = T)
sDMRs <- setdiff(sDMRs,one_kb_downstream,ignore.strand = T)
# Intergenic --------------------------------------------------------------
## Storing length of DMRs that are less that 1kb downstream of genes
cDMRs_overlap <- rbind(cDMRs_overlap,data.frame(Type = "Intergenic",Length = sum(width(cDMRs))))
sDMRs_overlap <- rbind(sDMRs_overlap,data.frame(Type = "Intergenic",Length = sum(width(sDMRs))))


####-####
###################################################################
#### Classifying DMRs based on overlap with other methods (10) ####
###################################################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs <- root_12h_CHH_ME_CE$`0.3`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE)

# Loading ATAC peaks -------------------------------------------------------
## Loading ATAC-seq peaks
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak")
ATAC_peaks <- keepSeqlevels(ATAC_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
ATAC_peaks <- GRanges(seqnames(ATAC_peaks),IRanges(start(ATAC_peaks),end(ATAC_peaks)))
ATAC_peaks <- unique(ATAC_peaks)
# Loading DNase peaks ------------------------------------------------------
## Loading DNase-seq peaks
DNase_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/DNase_q005_peaks.narrowPeak")
seqlevels(DNase_peaks) <- paste0("Chr",seqlevels(DNase_peaks))
DNase_peaks <- keepSeqlevels(DNase_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
DNase_peaks <- GRanges(seqnames(DNase_peaks),IRanges(start(DNase_peaks),end(DNase_peaks)))
DNase_peaks <- unique(DNase_peaks)
# Loading RNA data --------------------------------------------------------
#RNA_data <- read.table("/Volumes/PhD_one/NOME_Seq/features/RNA/Col0_ReadsPerGene.out.tab")
RNA_data <- read.table("/Volumes/shared/Shared210/Ryan/ecotypes_rna/RNA/counts/edited/Col0_1_DEseq.input")
#RNA_data <- RNA_data[-c(1:4),1:2]
# Loading GFF -------------------------------------------------------------
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Extracting gene sizes ---------------------------------------------------
exons_list_per_gene <- exonsBy(GFF,by="gene")
exonic_gene_sizes = list()
for (gene in names(exons_list_per_gene)) {
  exonic_gene_sizes[[gene]] = sum(width(reduce_ranges(exons_list_per_gene[[gene]])))
}
lengths = data.frame(id = names(exonic_gene_sizes), length = unlist(exonic_gene_sizes), row.names = NULL)
# Calculating TPM ---------------------------------------------------------
tpm <- RNA_data %>%
  rename(id = V1,count = V2) %>%
  inner_join(lengths) %>%
  mutate(tpm = count/length,tpm = tpm/sum(tpm)) %>%
  select(id,tpm)
# Classifying genes -------------------------------------------------------
GFF_genes <- genes_txdb(GFF)
gene_size_threshold <- 2000

classI_genes <- GFF_genes[overlapsAny(GFF_genes, DMRs) &  overlapsAny(GFF_genes, ATAC_peaks) & countOverlaps(GFF_genes,GFF_genes, maxgap = gene_size_threshold) == 1]
classII_genes <- GFF_genes[overlapsAny(GFF_genes, DMRs) &  !overlapsAny(GFF_genes, ATAC_peaks) & countOverlaps(GFF_genes,GFF_genes, maxgap = gene_size_threshold)== 1]
classIII_genes <- GFF_genes[!overlapsAny(GFF_genes, DMRs) &  overlapsAny(GFF_genes, ATAC_peaks) & countOverlaps(GFF_genes,GFF_genes, maxgap = gene_size_threshold)== 1]
classIV_genes <- GFF_genes[!overlapsAny(GFF_genes, DMRs) &  !overlapsAny(GFF_genes, ATAC_peaks) & countOverlaps(GFF_genes,GFF_genes, maxgap = gene_size_threshold)== 1]

classI_genes <- tpm[tpm$id %in% classI_genes$gene_id,] %>%
  mutate(class = "Common")

classII_genes <- tpm[tpm$id %in% classII_genes$gene_id,] %>%
  mutate(class = "iNOME")

classIII_genes <- tpm[tpm$id %in% classIII_genes$gene_id,] %>%
  mutate(class = "ATAC")

classIV_genes <- tpm[tpm$id %in% classIV_genes$gene_id,] %>%
  mutate(class = "None")

allClasses <- bind_rows(classI_genes,classII_genes,classIII_genes,classIV_genes)
####-####
###################################
#### Features around DMRs (11) ####
###################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs_root <- root_12h_CHH_ME_CE$`0.3`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE)
# Loading ATAC peaks -------------------------------------------------------
## Loading ATAC-seq peaks
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak")
ATAC_peaks <- keepSeqlevels(ATAC_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
ATAC_peaks <- GRanges(seqnames(ATAC_peaks),IRanges(start(ATAC_peaks),end(ATAC_peaks)))
ATAC_peaks <- unique(ATAC_peaks)
# Loading DNase peaks ------------------------------------------------------
## Loading DNase-seq peaks
DNase_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/DNase_q005_peaks.narrowPeak")
seqlevels(DNase_peaks) <- paste0("Chr",seqlevels(DNase_peaks))
DNase_peaks <- keepSeqlevels(DNase_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
DNase_peaks <- GRanges(seqnames(DNase_peaks),IRanges(start(DNase_peaks),end(DNase_peaks)))
DNase_peaks <- unique(DNase_peaks)
# Splitting DMRs ----------------------------------------------------------
## Extracting common DMRs
cDMRs <- intersect(DMRs_root,ATAC_peaks)
cDMRs <- intersect(cDMRs,DNase_peaks)
## Extracting unique DMRs
sDMRs <- setdiff(DMRs_root,ATAC_peaks)
sDMRs <- setdiff(sDMRs,DNase_peaks)
## Saving backup
cDMRs_backup <- cDMRs
sDMRs_backup <- sDMRs
# Loading Features --------------------------------------------------------
files <- dir("/Volumes/PhD_one/NOME_Seq/Col_features/bw/")
features <- vector("list",length(files))
for(i in 1:length(files)){
  features[[i]] <- import(paste0("/Volumes/PhD_one/NOME_Seq/Col_features/bw/",files[i]))
  features[[i]] <- keepSeqlevels(features[[i]],c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
}
names(features) <- gsub(pattern = "_col_10.bw",replacement = "",x = files)
names(features) <- gsub(pattern = "_root_edited.bw",replacement = "",x = names(features)) 
# Normalizing features around cDMR, sDMRs and all DMRs --------------------
normalised_DMRs <- vector("list",length(features))
names(normalised_DMRs) <- names(features)
pb <- txtProgressBar(0,length(features),style = 3)
for(i in 1:length(normalised_DMRs)){
  normalised_DMRs[[i]] <- normalizeToMatrix(features[[i]],DMRs_root,
                                            value_column = "score",
                                            mean_mode = "absolute",
                                            extend = 1000,w = 10,
                                            background = NA,smooth = F,
                                            target_ratio = 0.3)
  setTxtProgressBar(pb,i)
}

normalised_cDMRs <- vector("list",length(features))
names(normalised_cDMRs) <- names(features)
pb <- txtProgressBar(0,length(features),style = 3)
for(i in 1:length(normalised_cDMRs)){
  normalised_cDMRs[[i]] <- normalizeToMatrix(features[[i]],cDMRs,
                                             value_column = "score",
                                             mean_mode = "absolute",
                                             extend = 1000,w = 10,
                                             background = NA,smooth = F,
                                             target_ratio = 0.3)
  setTxtProgressBar(pb,i)
}

normalised_sDMRs <- vector("list",length(features))
names(normalised_sDMRs) <- names(features)
pb <- txtProgressBar(0,length(features),style = 3)
for(i in 1:length(normalised_sDMRs)){
  normalised_sDMRs[[i]] <- normalizeToMatrix(features[[i]],sDMRs,
                                             value_column = "score",
                                             mean_mode = "absolute",
                                             extend = 1000,w = 10,
                                             background = NA,smooth = F,
                                             target_ratio = 0.3)
  setTxtProgressBar(pb,i)
}
# Converting normalized matrix to data frame ------------------------------
normalised_DMRs_df <- data.frame(Position = numeric(),Score = numeric(),Feature = character(),Type = character())
for(i in 1:length(normalised_cDMRs)){
  normalised_DMRs_df <- rbind(normalised_DMRs_df,
                              data.frame("Position" = seq(10,2860,10),
                                         "Score" = colMeans(normalised_DMRs[[i]],na.rm = T),
                                         Feature = names(normalised_DMRs)[i],Type = "All Regions"))
}

normalised_cDMRs_df <- data.frame(Position = numeric(),Score = numeric(),Feature = character(),Type = character())
for(i in 1:length(normalised_cDMRs)){
  normalised_cDMRs_df <- rbind(normalised_cDMRs_df,
                               data.frame("Position" = seq(10,2860,10),
                                          "Score" = colMeans(normalised_cDMRs[[i]],na.rm = T),
                                          Feature = names(normalised_cDMRs)[i],Type = "Common Regions"))
}

normalised_sDMRs_df <- data.frame(Position = numeric(),Score = numeric(),Feature = character(),Type = character())
for(i in 1:length(normalised_sDMRs)){
  normalised_sDMRs_df <- rbind(normalised_sDMRs_df,
                               data.frame("Position" = seq(10,2860,10),
                                          "Score" = colMeans(normalised_sDMRs[[i]],na.rm = T),
                                          Feature = names(normalised_sDMRs)[i],Type = "Unique Regions"))
}

normalised_df <- rbind(normalised_DMRs_df,normalised_cDMRs_df,normalised_sDMRs_df)

####-####
##############################################
#### Features around intergenic DMRs (12) ####
##############################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs_root <- root_12h_CHH_ME_CE$`0.3`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE)
# Loading ATAC peaks -------------------------------------------------------
## Loading ATAC-seq peaks
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak")
ATAC_peaks <- keepSeqlevels(ATAC_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
ATAC_peaks <- GRanges(seqnames(ATAC_peaks),IRanges(start(ATAC_peaks),end(ATAC_peaks)))
ATAC_peaks <- unique(ATAC_peaks)
# Loading DNase peaks ------------------------------------------------------
## Loading DNase-seq peaks
DNase_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/DNase_q005_peaks.narrowPeak")
seqlevels(DNase_peaks) <- paste0("Chr",seqlevels(DNase_peaks))
DNase_peaks <- keepSeqlevels(DNase_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
DNase_peaks <- GRanges(seqnames(DNase_peaks),IRanges(start(DNase_peaks),end(DNase_peaks)))
DNase_peaks <- unique(DNase_peaks)
# Loading GFF -------------------------------------------------------------
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Loading GFF (Transposable elements) -------------------------------------
TE <- import("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
TE <- keepSeqlevels(TE,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
TE <- TE[TE$type=="transposable_element"]
# Loading Features --------------------------------------------------------
files <- dir("/Volumes/PhD_one/NOME_Seq/Col_features/bw/")
features <- vector("list",length(files))
for(i in 1:length(files)){
  features[[i]] <- import(paste0("/Volumes/PhD_one/NOME_Seq/Col_features/bw/",files[i]))
  features[[i]] <- keepSeqlevels(features[[i]],c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
}
names(features) <- gsub(pattern = "_col_10.bw",replacement = "",x = files)
names(features) <- gsub(pattern = "_root_edited.bw",replacement = "",x = names(features))

# Splitting DMRs ----------------------------------------------------------
## Extracting common DMRs
cDMRs <- intersect(DMRs_root,ATAC_peaks)
cDMRs <- intersect(cDMRs,DNase_peaks)
## Extracting unique DMRs
sDMRs <- setdiff(DMRs_root,ATAC_peaks)
sDMRs <- setdiff(sDMRs,DNase_peaks)
## Saving backup
cDMRs_backup <- cDMRs
sDMRs_backup <- sDMRs
# Extracting intergenic regions -------------------------------------------
DMRs_root <- setdiff(DMRs_root,genes_txdb(GFF),ignore.strand = T)
cDMRs <- setdiff(cDMRs,genes_txdb(GFF),ignore.strand = T)
sDMRs <- setdiff(sDMRs,genes_txdb(GFF),ignore.strand = T)

DMRs_root <- setdiff(DMRs_root,TE,ignore.strand = T)
cDMRs <- setdiff(cDMRs,TE,ignore.strand = T)
sDMRs <- setdiff(sDMRs,TE,ignore.strand = T)

DMRs_root <- setdiff(DMRs_root,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)
cDMRs <- setdiff(cDMRs,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)
sDMRs <- setdiff(sDMRs,promoters(GFF,upstream = 1000,downstream = 0),ignore.strand = T)

one_kb_downstream <- genes_txdb(GFF)
start(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"])
end(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"]) + 1000
end(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"])
start(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"]) - 1000

DMRs_root <- setdiff(DMRs_root,one_kb_downstream,ignore.strand = T)
cDMRs <- setdiff(cDMRs,one_kb_downstream,ignore.strand = T)
sDMRs <- setdiff(sDMRs,one_kb_downstream,ignore.strand = T)
# Normalizing features around cDMR, sDMRs and all DMRs --------------------
normalised_DMRs <- vector("list",length(features))
names(normalised_DMRs) <- names(features)
pb <- txtProgressBar(0,length(features),style = 3)
for(i in 1:length(normalised_DMRs)){
  normalised_DMRs[[i]] <- normalizeToMatrix(features[[i]],DMRs_root,
                                            value_column = "score",
                                            mean_mode = "absolute",
                                            extend = 1000,w = 10,
                                            background = NA,smooth = F,
                                            target_ratio = 0.3)
  setTxtProgressBar(pb,i)
}

normalised_cDMRs <- vector("list",length(features))
names(normalised_cDMRs) <- names(features)
pb <- txtProgressBar(0,length(features),style = 3)
for(i in 1:length(normalised_cDMRs)){
  normalised_cDMRs[[i]] <- normalizeToMatrix(features[[i]],cDMRs,
                                             value_column = "score",
                                             mean_mode = "absolute",
                                             extend = 1000,w = 10,
                                             background = NA,smooth = F,
                                             target_ratio = 0.3)
  setTxtProgressBar(pb,i)
}

normalised_sDMRs <- vector("list",length(features))
names(normalised_sDMRs) <- names(features)
pb <- txtProgressBar(0,length(features),style = 3)
for(i in 1:length(normalised_sDMRs)){
  normalised_sDMRs[[i]] <- normalizeToMatrix(features[[i]],sDMRs,
                                             value_column = "score",
                                             mean_mode = "absolute",
                                             extend = 1000,w = 10,
                                             background = NA,smooth = F,
                                             target_ratio = 0.3)
  setTxtProgressBar(pb,i)
}Ø
# Converting normalized matrix to data frame ------------------------------
normalised_DMRs_df <- data.frame(Position = numeric(),Score = numeric(),Feature = character(),Type = character())
for(i in 1:length(normalised_cDMRs)){
  normalised_DMRs_df <- rbind(normalised_DMRs_df,
                              data.frame("Position" = seq(10,2860,10),
                                         "Score" = colMeans(normalised_DMRs[[i]],na.rm = T),
                                         Feature = names(normalised_DMRs)[i],Type = "All Regions"))
}

normalised_cDMRs_df <- data.frame(Position = numeric(),Score = numeric(),Feature = character(),Type = character())
for(i in 1:length(normalised_cDMRs)){
  normalised_cDMRs_df <- rbind(normalised_cDMRs_df,
                               data.frame("Position" = seq(10,2860,10),
                                          "Score" = colMeans(normalised_cDMRs[[i]],na.rm = T),
                                          Feature = names(normalised_cDMRs)[i],Type = "Common Regions"))
}

normalised_sDMRs_df <- data.frame(Position = numeric(),Score = numeric(),Feature = character(),Type = character())
for(i in 1:length(normalised_sDMRs)){
  normalised_sDMRs_df <- rbind(normalised_sDMRs_df,
                               data.frame("Position" = seq(10,2860,10),
                                          "Score" = colMeans(normalised_sDMRs[[i]],na.rm = T),
                                          Feature = names(normalised_sDMRs)[i],Type = "Unique Regions"))
}

normalised_df <- rbind(normalised_DMRs_df,normalised_cDMRs_df,normalised_sDMRs_df)
####-####
##################################################
#### iNOME and other features around TSS (13,14) ####
##################################################
####-####
# Loading Features --------------------------------------------------------
files <- dir("/Volumes/PhD_one/NOME_Seq/Col_features/bw/")
features <- vector("list",length(files))
for(i in 1:length(files)){
  features[[i]] <- import(paste0("/Volumes/PhD_one/NOME_Seq/Col_features/bw/",files[i]))
  features[[i]] <- keepSeqlevels(features[[i]],c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
}
names(features) <- gsub(pattern = "_col_10.bw",replacement = "",x = files)
names(features) <- gsub(pattern = "_root_edited.bw",replacement = "",x = names(features))

# Loading RNA data --------------------------------------------------------
RNA_data <- read.table("/Volumes/PhD_one/NOME_Seq/features/RNA/Col0_ReadsPerGene.out.tab")
RNA_data <- RNA_data[-c(1:4),1:2]
# Loading GFF ------------------------------------------------------------- 
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Extracting gene sizes ---------------------------------------------------
exons_list_per_gene <- exonsBy(GFF,by="gene")
exonic_gene_sizes = list()
for (gene in names(exons_list_per_gene)) {
  exonic_gene_sizes[[gene]] = sum(width(reduce(exons_list_per_gene[[gene]])))
}
lengths = data.frame(id = names(exonic_gene_sizes), length = unlist(exonic_gene_sizes), row.names = NULL)
# Calculating TPM ---------------------------------------------------------
tpm <- data.frame(id = character(),TPM = numeric())
for(i in 1:nrow(RNA_data)){
  tpm <- rbind(tpm,data.frame(id = RNA_data[lengths[i,1]==RNA_data[,1],1],
                              TPM = RNA_data[lengths[i,1]==RNA_data[,1],2]/lengths[lengths[i,1]==RNA_data[,1],2]))
}
tpm$TPM <- tpm$TPM/(sum(tpm$TPM)/1000000)
# Loading TSS and ordering by expression ----------------------------------
tss <- TSS(GFF)
tss$TPM <- NA
pb <- txtProgressBar(0,nrow(tpm),style = 3)
for(i in 1:nrow(tpm)){
  if(length(tss[names(tss)==tpm[i,1]]$TPM)>0){
    tss[names(tss)==tpm[i,1]]$TPM <- tpm[i,2]
  }
  setTxtProgressBar(pb,i)
}
tss <- tss[order(tss$TPM,decreasing = T)]
tss$rank <- seq(1,length(tss),1)
genenames <- names(tss)
names(tss) <- tss$rank
tss <- dropSeqlevels(tss,c("ChrM","ChrC"),"coarse")
# Settings ----------------------------------------------------------------
ecotype <- "Col0" # Col0 , Cvi , Ler , C24 , WS2
movingwindow <- 100
# Loading data ------------------------------------------------------------
if(ecotype=="Col0"){
  files <- dir(path = paste0("~/work/Julia/cx_root_vs_leaf/root/"),pattern = "12h")
  cx_report <- lapply(paste0("~/work/Julia/cx_root_vs_leaf/root/",files),function(x){
    x <- readBismark2(x)
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  })
  names(cx_report) <- gsub(pattern = "_CX_report.txt",replacement = "",files)
} else{
  files <- dir(path = paste0("~/work/ecotypes/NOME/CX_reports/",ecotype,"/"))
  cx_report <- lapply(paste0("~/work/ecotypes/NOME/CX_reports/",ecotype,"/",files),function(x){
    x <- readBismark2(x)
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  })
  names(cx_report) <- gsub(pattern = "_R_fp_bismark_bt2_pe.deduplicated.CX_report.txt",replacement = "",files)
}
# Extracting GC -----------------------------------------------------------
if(ecotype=="Col0"){
  genome <- import(paste0("~/work/ecotypes/annotations/",ecotype,"/genome/",dir(paste0(path = "~/work/ecotypes/annotations/",ecotype,"/genome/"),pattern = "fasta")))[1:5]
  names(genome) <- gsub(pattern = " CHROMOSOME dumped from ADB: Feb/3/09 16:9; last updated: 2009-02-02",
                        replacement = "",names(genome))
  names(genome) <- gsub(pattern = " CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02",
                        replacement = "",names(genome))
} else{
  genome <- import(paste0("~/work/ecotypes/annotations/",ecotype,"/genome/",dir(paste0(path = "~/work/ecotypes/annotations/",ecotype,"/genome/"),pattern = "fasta")))[1:5]
}

GC <- lapply(cx_report,extractGC,genome = genome,contexts = "ALL")
# Calculating moving proportion -------------------------------------------
if(!file.exists(paste0("~/work/Julia/ecotypes_comparison/rdata/","Col0","_",movingwindow,".rdata"))){
  GC_prop <- vector("list",length(GC))
  names(GC_prop) <- names(GC)
  for(i in 1:length(GC)){
    message(names(GC)[i])
    pb <- txtProgressBar(0,length(GC[[i]]),style = 3)
    GC_prop[[i]] <-vector("list",length(GC[[i]]))
    names(GC_prop[[i]]) <- names(GC[[i]])
    for(j in 1:length(GC[[i]])){
      GC_prop[[i]][[j]] <- movingproportion(GC[[i]][[j]],100)
      GC_prop[[i]][[j]] <- GC_prop[[i]][[j]][seq(1,nrow(GC_prop[[i]][[j]]),10),]
      colnames(GC_prop[[i]][[j]]) <- c("chromosome","start","end","score")
      GC_prop[[i]][[j]] <- GRanges(seqnames = GC_prop[[i]][[j]]$chromosome,
                                   ranges = IRanges(GC_prop[[i]][[j]]$start,
                                                    (GC_prop[[i]][[j]]$end + 9)),
                                   score = GC_prop[[i]][[j]]$score)
      seqlengths(GC_prop[[i]][[j]]) <- end(getWholeChromosomes(GC_prop[[i]][[j]]))
      setTxtProgressBar(pb,j)
    }
    message("\n")
  }
  save(GC_prop,file = paste0("~/work/Julia/ecotypes_comparison/rdata/","Col0","_",movingwindow,".rdata"))
} else{
  load(paste0("~/work/Julia/ecotypes_comparison/rdata/","Col0","_",movingwindow,".rdata"))
}
# Computing moving iNOME signal (split) -----------------------------------
if(ecotype=="Col0"){
  control <- GC_prop[grep(pattern = "CE",x = names(GC_prop))][[1]]
  treatment <- GC_prop[grep(pattern = "ME",x = names(GC_prop))][[1]]
} else{
  control <- GC_prop[grep(pattern = "0h",x = names(GC_prop))][[1]]
  treatment <- GC_prop[grep(pattern = "12h",x = names(GC_prop))][[1]]
}

for(i in 1:length(treatment)){
  treatment[[i]]$score <- treatment[[i]]$score - control[[i]]$score
  treatment[[i]] <- treatment[[i]][!is.na(treatment[[i]]$score)]
  treatment[[i]] <- treatment[[i]][is.finite(treatment[[i]]$score)]
  treatment[[i]]$minus_score <- 1 - treatment[[i]]$score 
  treatment[[i]] <- treatment[[i]][!is.na(treatment[[i]]$minus_score)]
  treatment[[i]] <- treatment[[i]][is.finite(treatment[[i]]$minus_score)]
}
# Normalizing data around TSS  --------------------------------------------
mat_12h <- vector("list",length(treatment))
names(mat_12h) <- names(treatment)
pb <- txtProgressBar(0,length(treatment),style = 3)
for(i in 1:length(treatment)){
  mat_12h[[i]] <- normalizeToMatrix(treatment[[i]],tss[tss$TPM>0],value_column = "score",
                                          extend = c(500,1000), mean_mode = "absolute",
                                          w = 10)
  setTxtProgressBar(pb,i)
}

## Normalising data
iNOME_12h <- iNOME[grep("12h",names(iNOME))]
mat_12h_minus <- vector("list",length(iNOME_12h))
pb <- txtProgressBar(0,length(iNOME_12h[[1]]),style = 3)
for(i in 1:length(iNOME_12h[[1]])){
  mat_12h_minus[[i]] <- normalizeToMatrix(iNOME_12h[[1]][[i]],tss[tss$TPM>0],value_column = "score",
                                    extend = c(500,1000), mean_mode = "absolute",
                                    w = 10)
  setTxtProgressBar(pb,i)
}

mat_mNASE <- normalizeToMatrix(features[["mNASE"]],tss,value_column = "score",
                               extend = c(500,1000), mean_mode = "absolute",
                               w = 10)
## Converting to data frame
names(mat_12h_minus) <- c("GCG","GCHG","GCHH")
mat_12h_minus_df <- data.frame(Position = numeric(),Score = numeric(),Type = character())
for(i in 1:length(mat_12h_minus)){
  mat_12h_minus_df <- rbind(mat_12h_minus_df,
                              data.frame("Position" = seq(10,1500,10),
                                         "Score" = colMeans(mat_12h_minus[[i]],na.rm = T),
                                         Type = names(mat_12h_minus)[i]))
}

mat_mNASE_df <- data.frame(Position = seq(10,1500,10),Score = colMeans(mat_mNASE,na.rm = T),Type = "mNASE")
####-####
#####################################
#### Comparing Col0 and LER (15) ####
#####################################
####-####
# Loading RNA data --------------------------------------------------------
files <- dir("~/work/Julia/ecotypes_comparison/RNA/")
files <- c(files[grep("LER",files)],files[grep("Col",files)])
RNA_data <- lapply(paste0("~/work/Julia/ecotypes_comparison/RNA/",files),function(x){
  x <- read.table(x)[-c(1:4),1:2]
  colnames(x) <- c("Gene_ID","Count")
  x
})
names(RNA_data) <- gsub(pattern = "_ReadsPerGene.out.tab",replacement = "",x = files)
rm(files)
# Loading blast outputs ---------------------------------------------------
Ler_blast <- read.table("~/work/Julia/ecotypes_comparison/blast/Ler_blast.tab")[,c(1:2)]
Cvi_blast <- read.table("~/work/Julia/ecotypes_comparison/blast/Cvi_blast.tab")[,c(1:2)]
C24_blast <- read.table("~/work/Julia/ecotypes_comparison/blast/C24_blast.tab")[,c(1:2)]
WS2_blast <- read.table("~/work/Julia/ecotypes_comparison/blast/WS2_blast.tab")[,c(1:2)]
colnames(Ler_blast) <- c("Ler","Col0")
colnames(Cvi_blast) <- c("Cvi","Col0")
colnames(C24_blast) <- c("C24","Col0")
colnames(WS2_blast) <- c("WS2","Col0")
# Converting gene names ---------------------------------------------------
Ler_RNA <- RNA_data[grep("LER",names(RNA_data))]
Col_RNA <- RNA_data[grep("Col",names(RNA_data))]

for(i in 1:length(Ler_RNA)){
  Ler_RNA[[i]] <- merge(Ler_RNA[[i]],Ler_blast,by.x = "Gene_ID",by.y = "Ler")
  Ler_RNA[[i]]$Gene_ID <- Ler_RNA[[i]]$Col0
  Ler_RNA[[i]] <- Ler_RNA[[i]][,c(1:2)]
}

for(i in 1:length(Ler_RNA)){
  Col_RNA[[i]] <- merge(Col_RNA[[i]],Ler_blast,by.x = "Gene_ID",by.y = "Col0")
  Col_RNA[[i]] <- Col_RNA[[i]][,c(1:2)]
}
rm(RNA_data)
# Preparing data for DESEQ2 -----------------------------------------------
star_matrix <- data.frame(Gene_ID = Col_RNA[[1]]$Gene_ID,Col1 = Col_RNA[[1]]$Count)
star_matrix <- merge(star_matrix,Col_RNA[[2]],by = "Gene_ID")
star_matrix <- unique(star_matrix)
star_matrix <- merge(star_matrix,Ler_RNA[[1]],by = "Gene_ID")
star_matrix <- unique(star_matrix)
star_matrix <- merge(star_matrix,Ler_RNA[[2]],by = "Gene_ID")
star_matrix <- unique(star_matrix)

star_matrix$Gene_ID <- gsub("\\..*","",star_matrix$Gene_ID)

x <- split(star_matrix,star_matrix$Gene_ID)
for(i in 1:length(x)){
  if(nrow(x[[i]])>1){
    y <- sample(seq(1,nrow(x[[i]])),1)
    x[[i]] <- x[[i]][y,]
  }
}
star_matrix <- do.call("rbind",x)
star_matrix <- star_matrix[,-c(1)]
colnames(star_matrix) <- c("Col1","Col2","LER1","LER2")
# Calculating differential expression -------------------------------------
for(i in 1:ncol(star_matrix)){
  x <- data.frame(rownames(star_matrix),star_matrix[,i])
  write.table(x,file = paste0("~/work/Julia/ecotypes_comparison/",colnames(star_matrix[i]),".txt"),
              sep = "\t",quote = F,append = F,row.names = F,col.names = F)
}
sampleTable <- data.frame(sampleName = c("Col1","Col2","LER1","LER2"),fileName = c("Col1.txt","Col2.txt",
                                                                                   "LER1.txt","LER2.txt"),
                          condition = c("Col0","Col0","LER","LER"))

ddsHTSeq <- DESeq2::DESeqDataSetFromMatrix(countData = star_matrix,colData = sampleTable,
                                           design = ~ condition)
dds <- DESeq2::DESeq(ddsHTSeq)

res = DESeq2::results(dds, contrast=c("condition","LER","Col0"))
res = DESeq2::lfcShrink(dds,contrast = c("condition","LER","Col0"),res = res,type = 'normal')
res <- res[!is.na(res$padj),]
topGenes <- res[res$padj < 0.001,][c(2,6)]
topGenes$name <- gsub("\\..*","",rownames(topGenes))
# Loading GFFs ------------------------------------------------------------
## Loading Col0 GFF
Col0_GFF <- import("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
Col0_GFF <- keepSeqlevels(Col0_GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
Col0_GFF <- Col0_GFF[Col0_GFF$type=="gene"]
Col0_GFF <- GRanges(seqnames(Col0_GFF),IRanges(start(Col0_GFF),end(Col0_GFF)),strand = strand(Col0_GFF),name = Col0_GFF$Name)
## Loading Ler GFF
Ler_GFF <- import("~/work/ecotypes/annotations/LER/AT7213.gff3")
Ler_GFF <- keepSeqlevels(Ler_GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
Ler_GFF <- Ler_GFF[Ler_GFF$type=="gene"]
Ler_GFF <- GRanges(seqnames(Ler_GFF),IRanges(start(Ler_GFF),end(Ler_GFF)),strand = strand(Ler_GFF),name = Ler_GFF$ID)
## Loading Cvi GFF
Cvi_GFF <- import("~/work/ecotypes/annotations/Cvi/AT6911.gff3")
Cvi_GFF <- keepSeqlevels(Cvi_GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
Cvi_GFF <- Cvi_GFF[Cvi_GFF$type=="gene"]
Cvi_GFF <- GRanges(seqnames(Cvi_GFF),IRanges(start(Cvi_GFF),end(Cvi_GFF)),strand = strand(Cvi_GFF),name = Cvi_GFF$ID)
## Loading C24 GFF
C24_GFF <- import("~/work/ecotypes/annotations/C24/AT6906.gff3")
C24_GFF <- keepSeqlevels(C24_GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
C24_GFF <- C24_GFF[C24_GFF$type=="gene"]
C24_GFF <- GRanges(seqnames(C24_GFF),IRanges(start(C24_GFF),end(C24_GFF)),strand = strand(C24_GFF),name = C24_GFF$ID)
## Loading WS2 GFF
WS2_GFF <- import("~/work/ecotypes/annotations/WS2/WS27.gff3")
WS2_GFF <- keepSeqlevels(WS2_GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
WS2_GFF <- WS2_GFF[WS2_GFF$type=="gene"]
WS2_GFF <- GRanges(seqnames(WS2_GFF),IRanges(start(WS2_GFF),end(WS2_GFF)),strand = strand(WS2_GFF),name = WS2_GFF$ID)
## Preparing Ler GFF
Ler_GFF <- merge(as.data.frame(Ler_GFF),Ler_blast,by.x = "name",by.y = "Ler")
x <- split(Ler_GFF,Ler_GFF$Col0)
for(i in 1:length(x)){
  if(nrow(x[[i]])>1){
    y <- sample(seq(1,nrow(x[[i]])),1)
    x[[i]] <- x[[i]][y,]
  }
}
Ler_GFF <- do.call("rbind",x)
Ler_GFF <- GRanges(Ler_GFF$seqnames,IRanges(Ler_GFF$start,Ler_GFF$end),strand = Ler_GFF$strand,Name = Ler_GFF$Col0)
## Preparing Cvi GFF
Cvi_GFF <- merge(as.data.frame(Cvi_GFF),Cvi_blast,by.x = "name",by.y = "Cvi")
x <- split(Cvi_GFF,Cvi_GFF$Col0)
for(i in 1:length(x)){
  if(nrow(x[[i]])>1){
    y <- sample(seq(1,nrow(x[[i]])),1)
    x[[i]] <- x[[i]][y,]
  }
}
Cvi_GFF <- do.call("rbind",x)
Cvi_GFF <- GRanges(Cvi_GFF$seqnames,IRanges(Cvi_GFF$start,Cvi_GFF$end),strand = Cvi_GFF$strand,Name = Cvi_GFF$Col0)
## Preparing C24 GFF
C24_GFF <- merge(as.data.frame(C24_GFF),C24_blast,by.x = "name",by.y = "C24")
x <- split(C24_GFF,C24_GFF$Col0)
for(i in 1:length(x)){
  if(nrow(x[[i]])>1){
    y <- sample(seq(1,nrow(x[[i]])),1)
    x[[i]] <- x[[i]][y,]
  }
}
C24_GFF <- do.call("rbind",x)
C24_GFF <- GRanges(C24_GFF$seqnames,IRanges(C24_GFF$start,C24_GFF$end),strand = C24_GFF$strand,Name = C24_GFF$Col0)
## Preparing WS2 GFF
WS2_GFF <- merge(as.data.frame(WS2_GFF),WS2_blast,by.x = "name",by.y = "WS2")
x <- split(WS2_GFF,WS2_GFF$Col0)
for(i in 1:length(x)){
  if(nrow(x[[i]])>1){
    y <- sample(seq(1,nrow(x[[i]])),1)
    x[[i]] <- x[[i]][y,]
  }
}
WS2_GFF <- do.call("rbind",x)
WS2_GFF <- GRanges(WS2_GFF$seqnames,IRanges(WS2_GFF$start,WS2_GFF$end),strand = WS2_GFF$strand,Name = WS2_GFF$Col0)

width(Col0_GFF[Col0_GFF$name=="AT2G25450"])
width(Ler_GFF[Ler_GFF$Name=="AT2G25450"])
width(Cvi_GFF[Cvi_GFF$Name=="AT2G25450"])
width(C24_GFF[C24_GFF$Name=="AT2G25450"])
width(WS2_GFF[WS2_GFF$Name=="AT2G25450"])


x <- data.frame(Name = c("AT2G25450"),
                Col0 = width(Col0_GFF[Col0_GFF$name=="AT2G25450"]),
                Ler = width(Ler_GFF[Ler_GFF$Name=="AT2G25450"]),
                Cvi = width(Cvi_GFF[Cvi_GFF$Name=="AT2G25450"]),
                C24 = width(C24_GFF[C24_GFF$Name=="AT2G25450"]),
                WS2 = width(WS2_GFF[WS2_GFF$Name=="AT2G25450"]))


x <- data.frame(Ecotype = c("Col0","Ler","Cvi","C24","WS2"),Gene = "AT2G25450",
                Chromosome = c(seqnames(Col0_GFF[Col0_GFF$name=="AT2G25450"]),
                               seqnames(Ler_GFF[Ler_GFF$Name=="AT2G25450"]),
                               seqnames(Cvi_GFF[Cvi_GFF$Name=="AT2G25450"]),
                               seqnames(C24_GFF[C24_GFF$Name=="AT2G25450"]),
                               seqnames(WS2_GFF[WS2_GFF$Name=="AT2G25450"])),
                Start = c(start(Col0_GFF[Col0_GFF$name=="AT2G25450"]),
                          start(Ler_GFF[Ler_GFF$Name=="AT2G25450"]),
                          start(Cvi_GFF[Cvi_GFF$Name=="AT2G25450"]),
                          start(C24_GFF[C24_GFF$Name=="AT2G25450"]),
                          start(WS2_GFF[WS2_GFF$Name=="AT2G25450"])),
                End = c(end(Col0_GFF[Col0_GFF$name=="AT2G25450"]),
                        end(Ler_GFF[Ler_GFF$Name=="AT2G25450"]),
                        end(Cvi_GFF[Cvi_GFF$Name=="AT2G25450"]),
                        end(C24_GFF[C24_GFF$Name=="AT2G25450"]),
                        end(WS2_GFF[WS2_GFF$Name=="AT2G25450"])))

x$Start <- x$Start - 500
x$End <- x$End + 500
x$IGV <- paste0(x$Chromosome,":",x$Start,"-",x$End)
export(Col0_GFF,"~/work/Julia/GTF/Col0.gtf")
export(Ler_GFF,"~/work/Julia/GTF/Ler.gtf")
export(Cvi_GFF,"~/work/Julia/GTF/Cvi.gtf")
export(C24_GFF,"~/work/Julia/GTF/C24.gtf")
export(WS2_GFF,"~/work/Julia/GTF/WS2.gtf")
# Loading DMRs ------------------------------------------------------------
## Loading Col0 DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CpG_ME_CE_12h.RData")
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHG_ME_CE_12h.RData")
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
Col0_DMRs <- GRangesList("CG" = root_12h_CpG_ME_CE$`0.3`$`700`,"CHG" = root_12h_CHG_ME_CE$`0.3`$`700`,
                         "CHH" = root_12h_CHH_ME_CE$`0.3`$`700`)
rm(root_12h_CpG_ME_CE,root_12h_CHG_ME_CE,root_12h_CHH_ME_CE)
## Loading Ler DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/Ler_DMRs.RData")
# Extracting accessible genes ---------------------------------------------
## Overlapping Col0 DMRs
Col0_genes <- GRangesList()
for(i in 1:length(Col0_DMRs)){
  Col0_genes[[i]] <- Col0_GFF[overlapsAny(Col0_GFF,Col0_DMRs[[i]])]
}
names(Col0_genes) <- c("CG","CHG","CHH")
## Combining Col genes
x <- unlist(Col0_genes)
x <- split(x,x$name)
for(i in 1:length(x)){
  if(length(x[[i]])>1){
    y <- sample(seq(1,length(x[[i]])),1)
    x[[i]] <- x[[i]][y]
  }
}
Col0_genes <- unlist(x)
## Overlapping Ler DMRs
Ler_genes <- GRangesList()
for(i in 1:length(Ler_DMRs)){
  Ler_genes[[i]] <- Ler_GFF[overlapsAny(Ler_GFF,Ler_DMRs[[i]])]
}
names(Ler_genes) <- c("CG","CHG","CHH")
## Combining Ler genes
x <- unlist(Ler_genes)
x <- split(x,x$Name)
for(i in 1:length(x)){
  if(length(x[[i]])>1){
    y <- sample(seq(1,length(x[[i]])),1)
    x[[i]] <- x[[i]][y]
  }
}
Ler_genes <- unlist(x)

# Comparing genes between col0 and ler ------------------------------------
Common_genes <- Col0_genes[Col0_genes$name %in% Ler_genes$Name]
Col0_only <- Col0_genes[!(Col0_genes$name %in% Ler_genes$Name)]
Ler_only <- Ler_genes[!(Ler_genes$Name %in% Col0_genes$name)]
# Comparing accessible genes to differentially expressed genes ------------
Common_genes <- merge(as.data.frame(Common_genes),as.data.frame(topGenes),by.x = "name",by.y = "name")
Common_genes <- GRanges(Common_genes$seqnames,IRanges(Common_genes$start,Common_genes$end),
                        name = Common_genes$name,log2FoldChange = Common_genes$log2FoldChange,
                        padj = Common_genes$padj)
Common_genes$regulation <- NA
Common_genes[Common_genes$log2FoldChange>0]$regulation <- "upregulated"
Common_genes[Common_genes$log2FoldChange<0]$regulation <- "downregulated"

Col0_only <- merge(as.data.frame(Col0_only),as.data.frame(topGenes),by.x = "name",by.y = "name")
Col0_only <- GRanges(Col0_only$seqnames,IRanges(Col0_only$start,Col0_only$end),
                     name = Col0_only$name,log2FoldChange = Col0_only$log2FoldChange,
                     padj = Col0_only$padj)
Col0_only$regulation <- NA
Col0_only[Col0_only$log2FoldChange>0]$regulation <- "upregulated"
Col0_only[Col0_only$log2FoldChange<0]$regulation <- "downregulated"

Ler_only <- merge(as.data.frame(Ler_only),as.data.frame(topGenes),by.x = "Name",by.y = "name")
Ler_only <- GRanges(Ler_only$seqnames,IRanges(Ler_only$start,Ler_only$end),
                    name = Ler_only$Name,log2FoldChange = Ler_only$log2FoldChange,
                    padj = Ler_only$padj)
Ler_only$regulation <- NA
Ler_only[Ler_only$log2FoldChange>0]$regulation <- "upregulated"
Ler_only[Ler_only$log2FoldChange<0]$regulation <- "downregulated"

# Building dataframe to store regulation data -----------------------------
regulation_df <- data.frame(Regulation = character(),Value = numeric(),Type = character())
regulation_df <- rbind(regulation_df,data.frame(Regulation = c("upregulated","downregulated"),
                                                Value = c(length(Common_genes[Common_genes$regulation=="upregulated"]),
                                                          length(Common_genes[Common_genes$regulation=="downregulated"])),
                                                Type = "Common"))
regulation_df <- rbind(regulation_df,data.frame(Regulation = c("upregulated","downregulated"),
                                                Value = c(length(Col0_only[Col0_only$regulation=="upregulated"]),
                                                          length(Col0_only[Col0_only$regulation=="downregulated"])),
                                                Type = "Col0 Only"))
regulation_df <- rbind(regulation_df,data.frame(Regulation = c("upregulated","downregulated"),
                                                Value = c(length(Ler_only[Ler_only$regulation=="upregulated"]),
                                                          length(Ler_only[Ler_only$regulation=="downregulated"])),
                                                Type = "Ler Only"))

x <- data.frame(c(NA,NA))
x <- cbind(x,"Common" = regulation_df[regulation_df$Type=="Common",2],"Col0 Only" = regulation_df[regulation_df$Type=="Col0 Only",2],
           "Ler Only" = regulation_df[regulation_df$Type=="Ler Only",2])[,2:4]
y <- data.frame(x[,1],x[,2])
fisher_Com_Col0 <- data.frame(fisher.test(y)$p.value)
y <- data.frame(x[,1],x[,3])
fisher_Com_LER <- data.frame(fisher.test(y)$p.value)
y <- data.frame(x[,2],x[,3])
fisher_Col0_LER <- data.frame(fisher.test(y)$p.value)
fisher_pvalue <- rbind("Common vs Col0" = fisher_Com_Col0,"Common vs LER" = fisher_Com_LER,"Col0 vs LER" = fisher_Col0_LER)
fisher_pvalue$type <- rownames(fisher_pvalue)
colnames(fisher_pvalue) <- c("Pvalue","Type")

ggplot(fisher_pvalue,aes(x = Type,y = Pvalue)) +
  geom_point() +
  geom_hline(yintercept = 0.05,linetype = "dotted",col = "red") +
  theme_cowplot() #(edit this!!!!!!!!!!!!!!!!!!)
####-####
###############################################
#### Comparing iNOME-seq to MNase-seq (16) ####
###############################################
####-####
# Loading RNA data --------------------------------------------------------
RNA_data <- read.table("/Volumes/PhD_one/NOME_Seq/features/RNA/Col0_ReadsPerGene.out.tab")
RNA_data <- RNA_data[-c(1:4),1:2]
# Loading GFF ------------------------------------------------------------- 
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Extracting gene sizes ---------------------------------------------------
exons_list_per_gene <- exonsBy(GFF,by="gene")
exonic_gene_sizes = list()
for (gene in names(exons_list_per_gene)) {
  exonic_gene_sizes[[gene]] = sum(width(reduce(exons_list_per_gene[[gene]])))
}
lengths = data.frame(id = names(exonic_gene_sizes), length = unlist(exonic_gene_sizes), row.names = NULL)
# Calculating TPM ---------------------------------------------------------
tpm <- data.frame(id = character(),TPM = numeric())
for(i in 1:nrow(RNA_data)){
  tpm <- rbind(tpm,data.frame(id = RNA_data[lengths[i,1]==RNA_data[,1],1],
                              TPM = RNA_data[lengths[i,1]==RNA_data[,1],2]/lengths[lengths[i,1]==RNA_data[,1],2]))
}
tpm$TPM <- tpm$TPM/(sum(tpm$TPM)/1000000)
# Loading TSS and ordering by expression ----------------------------------
tss <- genes_txdb(GFF)
end(tss[strand(tss)=="+"]) <- start(tss[strand(tss)=="+"])
start(tss[strand(tss)=="-"]) <- end(tss[strand(tss)=="-"])
tss$TPM <- NA
pb <- txtProgressBar(0,nrow(tpm),style = 3)
for(i in 1:nrow(tpm)){
  if(length(tss[names(tss)==tpm[i,1]]$TPM)>0){
    tss[names(tss)==tpm[i,1]]$TPM <- tpm[i,2]
  }
  setTxtProgressBar(pb,i)
}
tss <- tss[order(tss$TPM,decreasing = T)]
tss$rank <- seq(1,length(tss),1)
genenames <- names(tss)
names(tss) <- tss$rank
tss <- dropSeqlevels(tss,c("ChrM","ChrC"),"coarse")
# Preparing data for plotting ---------------------------------------------
quartiles <- 3
tss$quartile <- ntile(tss$TPM,quartiles)

mNASE_bw <- vector("list",quartiles)
iNOME_bw <- vector("list",quartiles)
for(i in 1:quartiles){
  mNASE_bw[[i]] <- regionPlot(bamFile = "~/work/Julia/ecotypes_comparison/mnase/mNASE_col_50.bw",testRanges = tss[tss$quartile==i],format = "bigwig",paired = T,distanceUp = 500,distanceDown = 1000)
  iNOME_bw[[i]] <- regionPlot(bamFile = "~/work/Julia/ecotypes_comparison/moving_average/Col0_CHH_100.bw",testRanges = tss[tss$quartile==i],format = "bigwig",paired = T,distanceUp = 500,distanceDown = 1000)
}

MNase_df <- vector("list",quartiles)
iNOME_df <- vector("list",quartiles)
for(i in 1:quartiles){
  MNase_df[[i]] <- assays(mNASE_bw[[i]])[[1]]
  iNOME_df[[i]] <- assays(iNOME_bw[[i]])[[1]]
}

for(i in 1:quartiles){
  MNase_df[[i]] <- data.frame(Position = seq(1,1501,1),Score = colMeans(MNase_df[[i]]),Type = "MNase",quartile = i)
  iNOME_df[[i]] <- data.frame(Position = seq(1,1501,1),Score = (1 - colMeans(iNOME_df[[i]]))-0.724,Type = "iNOME",quartile = i)
}
MNase_df <- do.call("rbind",MNase_df)
iNOME_df <- do.call("rbind",iNOME_df)
if(length(unique(MNase_df$quartile))==3){
  MNase_df[MNase_df$quartile==1,4] <- "Low"
  MNase_df[MNase_df$quartile==2,4] <- "Medium"
  MNase_df[MNase_df$quartile==3,4] <- "High"
  
  iNOME_df[iNOME_df$quartile==1,4] <- "Low"
  iNOME_df[iNOME_df$quartile==2,4] <- "Medium"
  iNOME_df[iNOME_df$quartile==3,4] <- "High"
  
  MNase_df$quartile <- factor(MNase_df$quartile,levels = c("High","Medium","Low"))
  iNOME_df$quartile <- factor(iNOME_df$quartile,levels = c("High","Medium","Low"))
}

iNOME_MNase_comparison <- rbind(MNase_df,iNOME_df)

# Calculating correlations ------------------------------------------------
results <- data.frame(Comparison = character(),p.value = numeric(),R = numeric())
for(i in c("Low","Medium","High")){
  x <- cor.test(iNOME_MNase_comparison[iNOME_MNase_comparison$quartile==i & iNOME_MNase_comparison$Type=="MNase",2],
                iNOME_MNase_comparison[iNOME_MNase_comparison$quartile==i & iNOME_MNase_comparison$Type=="iNOME",2],
                method = "pearson")
  results <- rbind(results,data.frame(Comparison = i,p.value = x$p.value,R = x$estimate))
}
####-####
###########################################################
#### Comparing accessibility to expression in ecotypes ####
###########################################################
####-####
# Loading data ------------------------------------------------------------
## Loading annotations
ecotypes <- c("Col0","Ler","Cvi","C24","WS2")
GFF <- list()
pb <- txtProgressBar(0,length(ecotypes),style = 3)
for(i in 1:length(ecotypes)){
  GFF[[i]] <- import.gff(paste0("~/work/ecotypes/annotations/",ecotypes[i],"/",dir(path = paste0("~/work/ecotypes/annotations/",ecotypes[i]),pattern = "gff")))
  setTxtProgressBar(pb,i)
}
names(GFF) <- ecotypes

GFF <- lapply(GFF, function(x){
  x <- x[x$type=="gene"]
  x <- GRanges(seqnames = seqnames(x),ranges = IRanges(start(x),end(x)),strand = strand(x),ID = x$ID)
  x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  x
})
## Loading genomes
genomes <- list()
for(i in 1:length(ecotypes)){
  files <- dir(paste0("~/work/Julia/fasta/",ecotypes[[i]]))
  genomes[[i]] <- import(paste0("~/work/Julia/fasta/",ecotypes[[i]],"/",files))
}
names(genomes) <- ecotypes
## Loading cx_reports
cx_reports <- list()
for(i in 1:length(ecotypes)){
  if(ecotypes[[i]]=="Col0"){
    files <- dir(path = paste0("~/work/Julia/cx_root_vs_leaf/root/"),pattern = "12h")
    cx_reports[[i]] <- lapply(paste0("~/work/Julia/cx_root_vs_leaf/root/",files),function(x){
      x <- readBismark2(x)
      x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
      x <- extractGC(x,genomes[[i]],"ALL")
      x
    })
    names(cx_reports[[i]]) <- gsub(pattern = "_CX_report.txt",replacement = "",files)
  } else{
    files <- dir(path = paste0("~/work/ecotypes/NOME/CX_reports/",ecotypes[[i]],"/"))
    cx_reports[[i]] <- lapply(paste0("~/work/ecotypes/NOME/CX_reports/",ecotypes[[i]],"/",files),function(x){
      x <- readBismark2(x)
      x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
      x <- extractGC(x,genomes[[i]],"ALL")
      x
    })
    names(cx_reports[[i]]) <- gsub(pattern = "_R_fp_bismark_bt2_pe.deduplicated.CX_report.txt",replacement = "",files)
    names(cx_reports[[i]]) <- gsub(pattern = "MCViP_",replacement = "",x = names(cx_reports[[i]]))
  }
}

cx_reports <- lapply(cx_reports,lapply,function(x){
  x <- unlist(x)
  x$context <- "CHH"
  x
})
names(cx_reports) <- ecotypes

# Computing average Accessibility for genes -------------------------------
average_Access <- list()
for(i in 1:length(cx_reports)){
  average_Access[[i]] <- list()
  for(j in  1:length(cx_reports[[i]])){
    average_Access[[i]][[j]] <- analyseReadsInsideRegionsForCondition(GFF[[i]],cx_reports[[i]][[j]],"CHH",cores = 4)
  }
  names(average_Access[[i]]) <- names(cx_reports[[i]])
}
names(average_Access) <- names(cx_reports)

for(i in 1:length(average_Access)){
  average_Access[[i]][[2]]$proportionCHH <- average_Access[[i]][[2]]$proportionCHH - average_Access[[i]][[1]]$proportionCHH
  average_Access[[i]] <- GRanges(seqnames(average_Access[[i]][[2]]),IRanges(start(average_Access[[i]][[2]]),end(average_Access[[i]][[2]])),
                                 strand = strand(average_Access[[i]][[2]]),ID = average_Access[[i]][[2]]$ID,proportion = average_Access[[i]][[2]]$proportionCHH)
  average_Access[[i]] <- as.data.frame(average_Access[[i]])
} 
rm("cx_reports","files","genomes","GFF","i","j","pb")

# Loading RNA data --------------------------------------------------------
ecotypes <- c("Col0","Ler","Cvi","C24","WS2")
RNA <- list()
pb <- txtProgressBar(0,length(ecotypes),style = 3)
for(i in 1:length(ecotypes)){
  RNA[[i]] <- read.table(paste0("~/work/Julia/ecotypes_comparison/Accessibility_comp/RNA/",dir(path = paste0("~/work/Julia/ecotypes_comparison/Accessibility_comp/RNA/"),pattern = ecotypes[i])))[-c(1:4),1:2]
  setTxtProgressBar(pb,i)
}
names(RNA) <- ecotypes
# Loading genes ----------------------------------------------------------
## Loading GFF
ecotypes <- c("Col0","Ler","Cvi","C24","WS2")
genes <- list()
for(i in 1:length(ecotypes)){
  if(ecotypes[[i]]=="Col0"){
    genes[[i]] <- makeTxDbFromGFF(paste0("~/work/ecotypes/annotations/",ecotypes[i],"/",dir(path = paste0("~/work/ecotypes/annotations/",ecotypes[i]),pattern = "gff")))
  } else{
    genes[[i]] <- import(paste0("~/work/ecotypes/annotations/",ecotypes[i],"/",dir(path = paste0("~/work/ecotypes/annotations/",ecotypes[i]),pattern = "gff")))
  }
}
names(genes) <- ecotypes
## Removing irrelevent data
exons_list_per_gene <- list()
pb <- txtProgressBar(0,length(ecotypes),style = 3,width = 40)
for(i in 1:length(ecotypes)){
  if(ecotypes[i]=="Col0"){
    exons_list_per_gene[[i]] <- keepSeqlevels(genes[[i]],c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
    exons_list_per_gene[[i]] <- exonsBy(genes[[i]],"gene")
  } else{
    exons_list_per_gene[[i]] <- keepSeqlevels(genes[[i]],c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
    exons_list_per_gene[[i]] <- exons_list_per_gene[[i]][exons_list_per_gene[[i]]$type=="exon"]
    exons_list_per_gene[[i]] <- GRanges(seqnames(exons_list_per_gene[[i]]),IRanges(start(exons_list_per_gene[[i]]),end(exons_list_per_gene[[i]])),
                                        strand(exons_list_per_gene[[i]]),gene_id = exons_list_per_gene[[i]]$gene_id)
    exons_list_per_gene[[i]] <- split(exons_list_per_gene[[i]],exons_list_per_gene[[i]]$gene_id)
  }
  setTxtProgressBar(pb,i)
}
names(exons_list_per_gene) <- ecotypes
x <- GRangesList() 
pb <- txtProgressBar(0,length(exons_list_per_gene[[1]]),style = 3,width = 40)
for(i in 1:length(exons_list_per_gene[[1]])){
  x[[i]] <- exons_list_per_gene[[1]][[i]][,-c(1,2)]
  x[[i]]$gene_name <- names(exons_list_per_gene[[1]][i])
  setTxtProgressBar(pb,i)
}
names(x) <- names(exons_list_per_gene[[1]])

exons_list_per_gene[[1]] <- x
# Extracting gene sizes ---------------------------------------------------
lengths <- list()
for(i in 1:length(exons_list_per_gene)){
  message(ecotypes[i])
  pb <- txtProgressBar(0,length(exons_list_per_gene[[i]]),style = 3,width = 40)
  exonic_gene_sizes <- c()
  for(j in 1:length(exons_list_per_gene[[i]])){
    exonic_gene_sizes <- c(exonic_gene_sizes,sum(width(reduce(exons_list_per_gene[[i]][[j]]))))
    setTxtProgressBar(pb,j)
  }
  lengths[[i]] <- data.frame(id = names(exons_list_per_gene[[i]]),length = exonic_gene_sizes,row.names = NULL)
  message("\n")
}

# Calculating TPM ---------------------------------------------------------
tpm <- list()
for(i in 1:length(RNA)){
  message(ecotypes[[i]])
  pb <- txtProgressBar(0,nrow(RNA[[i]]),style = 3,width = 40)
  tpm[[i]] <- data.frame(id = character(),TPM = numeric())
  for(j in 1:nrow(RNA[[i]])){
    tpm[[i]] <- rbind(tpm[[i]],data.frame(id = RNA[[i]][lengths[[i]][j,1]==RNA[[i]][,1],1],
                                          TPM = RNA[[i]][lengths[[i]][j,1]==RNA[[i]][,1],2]/lengths[[i]][lengths[[i]][j,1]==RNA[[i]][,1],2]))
    setTxtProgressBar(pb,j)
  }
  tpm[[i]]$TPM <- tpm[[i]]$TPM/(sum(tpm[[i]]$TPM,na.rm = T)/1000000)
  message("\n")
}
names(tpm) <- ecotypes

# Merging accessibility and expression ------------------------------------
for(i in 1:length(average_Access)){
  average_Access[[i]] <- merge(average_Access[[i]],tpm[[i]],by.x = "ID",by.y = "id")
  average_Access[[i]] <- makeGRangesFromDataFrame(average_Access[[i]],keep.extra.columns = T)
  average_Access[[i]] <- average_Access[[i]][!is.na(average_Access[[i]]$TPM)]
  average_Access[[i]] <- average_Access[[i]][!is.na(average_Access[[i]]$proportion)]
}

####-####
#########################################
#### Calculating polII pausing index ####
#########################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs <- root_12h_CHH_ME_CE$`0.3`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE)

# Loading ATAC peaks -------------------------------------------------------
## Loading ATAC-seq peaks
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak")
ATAC_peaks <- keepSeqlevels(ATAC_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
ATAC_peaks <- GRanges(seqnames(ATAC_peaks),IRanges(start(ATAC_peaks),end(ATAC_peaks)))
ATAC_peaks <- unique(ATAC_peaks)
# Loading DNase peaks ------------------------------------------------------
## Loading DNase-seq peaks
DNase_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/DNase_q005_peaks.narrowPeak")
seqlevels(DNase_peaks) <- paste0("Chr",seqlevels(DNase_peaks))
DNase_peaks <- keepSeqlevels(DNase_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
DNase_peaks <- GRanges(seqnames(DNase_peaks),IRanges(start(DNase_peaks),end(DNase_peaks)))
DNase_peaks <- unique(DNase_peaks)
# Loading RNA data --------------------------------------------------------
RNA_data <- read.table("/Volumes/PhD_one/NOME_Seq/features/RNA/Col0_ReadsPerGene.out.tab")
RNA_data <- RNA_data[-c(1:4),1:2]
# Loading GFF -------------------------------------------------------------
## Loading GFF
GFF <- makeTxDbFromGFF("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff")
## Removing irrelevent data
GFF <- keepSeqlevels(GFF,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
# Extracting gene sizes ---------------------------------------------------
exons_list_per_gene <- exonsBy(GFF,by="gene")
exonic_gene_sizes = list()
for (gene in names(exons_list_per_gene)) {
  exonic_gene_sizes[[gene]] = sum(width(reduce(exons_list_per_gene[[gene]])))
}
lengths = data.frame(id = names(exonic_gene_sizes), length = unlist(exonic_gene_sizes), row.names = NULL)
# Calculating TPM ---------------------------------------------------------
tpm <- data.frame(id = character(),TPM = numeric())
for(i in 1:nrow(RNA_data)){
  tpm <- rbind(tpm,data.frame(id = RNA_data[lengths[i,1]==RNA_data[,1],1],
                              TPM = RNA_data[lengths[i,1]==RNA_data[,1],2]/
                                lengths[lengths[i,1]==RNA_data[,1],2]))
}
tpm$TPM <- tpm$TPM/sum(tpm$TPM)
# Classifying genes -------------------------------------------------------
GFF_genes <- genes(GFF)
gene_size_threshold <- 2000

classI_genes <- GFF_genes[overlapsAny(GFF_genes, DMRs) &  overlapsAny(GFF_genes, ATAC_peaks) & countOverlaps(GFF_genes,GFF_genes, maxgap = gene_size_threshold) == 1]
classII_genes <- GFF_genes[overlapsAny(GFF_genes, DMRs) &  !overlapsAny(GFF_genes, ATAC_peaks) & countOverlaps(GFF_genes,GFF_genes, maxgap = gene_size_threshold)== 1]
classIII_genes <- GFF_genes[!overlapsAny(GFF_genes, DMRs) &  overlapsAny(GFF_genes, ATAC_peaks) & countOverlaps(GFF_genes,GFF_genes, maxgap = gene_size_threshold)== 1]
classIV_genes <- GFF_genes[!overlapsAny(GFF_genes, DMRs) &  !overlapsAny(GFF_genes, ATAC_peaks) & countOverlaps(GFF_genes,GFF_genes, maxgap = gene_size_threshold)== 1]

classI_genes$TPM <- NA
pb <- txtProgressBar(0,length(classI_genes),style = 3)
for(i in 1:length(classI_genes)){
  classI_genes$TPM[i] <- tpm[classI_genes$gene_id[i]==tpm$id,2]
  setTxtProgressBar(pb,i)
}

classII_genes$TPM <- NA
pb <- txtProgressBar(0,length(classII_genes),style = 3)
for(i in 1:length(classII_genes)){
  classII_genes$TPM[i] <- tpm[classII_genes$gene_id[i]==tpm$id,2]
  setTxtProgressBar(pb,i)
}

classIII_genes$TPM <- NA
pb <- txtProgressBar(0,length(classIII_genes),style = 3)
for(i in 1:length(classIII_genes)){
  classIII_genes$TPM[i] <- tpm[classIII_genes$gene_id[i]==tpm$id,2]
  setTxtProgressBar(pb,i)
}

classIV_genes$TPM <- NA
pb <- txtProgressBar(0,length(classIV_genes),style = 3)
for(i in 1:length(classIV_genes)){
  classIV_genes$TPM[i] <- tpm[classIV_genes$gene_id[i]==tpm$id,2]
  setTxtProgressBar(pb,i)
}

all_classes <- data.frame(class = "class 1",TPM = classI_genes$TPM )
all_classes <- rbind(all_classes,data.frame(class = "class 2",TPM = classII_genes$TPM ))
all_classes <- rbind(all_classes,data.frame(class = "class 3",TPM = classIII_genes$TPM ))
all_classes <- rbind(all_classes,data.frame(class = "class 4",TPM = classIV_genes$TPM ))

####-####
###############################################
#### Comparing genes in different ecotypes ####
###############################################
####-####
# Settings ----------------------------------------------------------------
ecotype <- "Col0" # Col0 , Cvi , Ler , C24 , WS2
movingwindow <- 100
# Loading data ------------------------------------------------------------
if(ecotype=="Col0"){
  files <- dir(path = paste0("~/work/Julia/cx_root_vs_leaf/root/"),pattern = "12h")
  cx_report <- lapply(paste0("~/work/Julia/cx_root_vs_leaf/root/",files),function(x){
    x <- readBismark2(x)
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  })
  names(cx_report) <- gsub(pattern = "_CX_report.txt",replacement = "",files)
} else{
  files <- dir(path = paste0("~/work/ecotypes/NOME/CX_reports/",ecotype,"/"))
  cx_report <- lapply(paste0("~/work/ecotypes/NOME/CX_reports/",ecotype,"/",files),function(x){
    x <- readBismark2(x)
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  })
  names(cx_report) <- gsub(pattern = "_R_fp_bismark_bt2_pe.deduplicated.CX_report.txt",replacement = "",files)
}
# Extracting GC -----------------------------------------------------------
if(ecotype=="Col0"){
  genome <- import(paste0("~/work/ecotypes/annotations/",ecotype,"/genome/",dir(paste0(path = "~/work/ecotypes/annotations/",ecotype,"/genome/"),pattern = "fasta")))[1:5]
  names(genome) <- gsub(pattern = " CHROMOSOME dumped from ADB: Feb/3/09 16:9; last updated: 2009-02-02",
                        replacement = "",names(genome))
  names(genome) <- gsub(pattern = " CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02",
                        replacement = "",names(genome))
} else{
  genome <- import(paste0("~/work/ecotypes/annotations/",ecotype,"/genome/",dir(paste0(path = "~/work/ecotypes/annotations/",ecotype,"/genome/"),pattern = "fasta")))[1:5]
}

GC <- lapply(cx_report,extractGC,genome = genome,contexts = "ALL")
####-####
# Loading blast outputs ---------------------------------------------------
blast <- read.table(paste0("~/work/Julia/ecotypes_comparison/blast/",ecotype,"_blast.tab"))[,c(1:2)]
colnames(blast) <- c("Ecotype","Col0")
blast <- blast[!duplicated(blast),]
blast <- split(blast,blast$Col0)
pb <- txtProgressBar(0,length(blast),style = 3,width = 50)
for(i in 1:length(blast)){
  if(nrow(blast[[i]])>1){
    blast[[i]] <- blast[[i]][1,]
  }
  setTxtProgressBar(pb,i)
}
blast <- do.call("rbind",blast)
# Editing GFF (Not required for Col0) -------------------------------------
GFF <- import.gff(paste0("~/work/ecotypes/annotations/",ecotype,"/",dir(path = paste0("~/work/ecotypes/annotations/",ecotype),pattern = "gff")))
GFF_gene <- GFF[GFF$type=="gene"]
GFF_gene <- GFF_gene[GFF_gene$Name!="None"]
GFF_gene <- as.data.frame(GFF_gene)
GFF_gene <- merge(GFF_gene,blast,by.x = "ID",by.y = "Ecotype")
GFF_gene$ID <- GFF_gene$Col0
GFF_gene$Name <- GFF_gene$Col0
GFF_gene <- GFF_gene[,1:35]

GFF_features <- GFF[GFF$type!="gene"]
GFF_features <- GFF_features[GFF_features$Name!="None"]
GFF_features <- as.data.frame(GFF_features)
GFF_features <- merge(GFF_features,blast,by.x = "gene_id",by.y = "Ecotype")
GFF_features$gene_id <- GFF_features$Col0
GFF_features$Name <- GFF_features$Col0
GFF_features <- GFF_features[,1:35]

GFF <- rbind(GFF_gene,GFF_features)
GFF <- makeGRangesFromDataFrame(GFF,keep.extra.columns = T)

export.gff3(GFF,paste0("~/work/Julia/ecotypes_comparison/GFFs//",ecotype,"_edited.gff3")) # Remember to sort this file after export
####-####
# Computing moving average (split) ----------------------------------------
if(!file.exists(paste0("~/work/Julia/ecotypes_comparison/rdata/",ecotype,"_",movingwindow,".rdata"))){
  GC_prop <- vector("list",length(GC))
  names(GC_prop) <- names(GC)
  for(i in 1:length(GC)){
    message(names(GC)[i])
    pb <- txtProgressBar(0,length(GC[[i]]),style = 3)
    GC_prop[[i]] <-vector("list",length(GC[[i]]))
    names(GC_prop[[i]]) <- names(GC[[i]])
    for(j in 1:length(GC[[i]])){
      GC_prop[[i]][[j]] <- movingproportion(GC[[i]][[j]],movingwindow)
      GC_prop[[i]][[j]] <- GC_prop[[i]][[j]][seq(1,nrow(GC_prop[[i]][[j]]),10),]
      colnames(GC_prop[[i]][[j]]) <- c("chromosome","start","end","score")
      GC_prop[[i]][[j]] <- GRanges(seqnames = GC_prop[[i]][[j]]$chromosome,
                                   ranges = IRanges(GC_prop[[i]][[j]]$start,
                                                    (GC_prop[[i]][[j]]$end + 9)),
                                   score = GC_prop[[i]][[j]]$score)
      seqlengths(GC_prop[[i]][[j]]) <- end(getWholeChromosomes(GC_prop[[i]][[j]]))
      setTxtProgressBar(pb,j)
    }
    message("\n")
  }
  save(GC_prop,file = paste0("~/work/Julia/ecotypes_comparison/rdata/",ecotype,"_",movingwindow,".rdata"))
} else{
  load(paste0("~/work/Julia/ecotypes_comparison/rdata/",ecotype,"_",movingwindow,".rdata"))
}

# Computing moving iNOME signal (split) -----------------------------------
if(ecotype=="Col0"){
  control <- GC_prop[grep(pattern = "CE",x = names(GC_prop))][[1]]
  treatment <- GC_prop[grep(pattern = "ME",x = names(GC_prop))][[1]]
} else{
  control <- GC_prop[grep(pattern = "0h",x = names(GC_prop))][[1]]
  treatment <- GC_prop[grep(pattern = "12h",x = names(GC_prop))][[1]]
}

for(i in 1:length(treatment)){
  treatment[[i]]$score <- treatment[[i]]$score - control[[i]]$score
  treatment[[i]] <- treatment[[i]][!is.na(treatment[[i]]$score)]
  treatment[[i]] <- treatment[[i]][is.finite(treatment[[i]]$score)]
}
# Exporting as bigwig (split) ---------------------------------------------
export(treatment[[1]],paste0("~/work/Julia/ecotypes_comparison/moving_average/",ecotype,"_",names(treatment)[1],"_",movingwindow,".bw"),"bw")
export(treatment[[2]],paste0("~/work/Julia/ecotypes_comparison/moving_average/",ecotype,"_",names(treatment)[2],"_",movingwindow,".bw"),"bw")
export(treatment[[3]],paste0("~/work/Julia/ecotypes_comparison/moving_average/",ecotype,"_",names(treatment)[3],"_",movingwindow,".bw"),"bw")
####-####
# Computing moving average ------------------------------------------------
GC_all <- GC
for(i in 1:length(GC_all)){
  GC_all[[i]] <- unlist(GC_all[[i]])
  GC_all[[i]]$context <- "CHH"
}

if(!file.exists(paste0("~/work/Julia/ecotypes_comparison/rdata/",ecotype,"_",movingwindow,"_all.rdata"))){
  pb <- txtProgressBar(0,length(GC_all),style = 3)
  GC_prop_all <- vector("list",length(GC_all))
  names(GC_prop_all) <- names(GC_all)
  for(i in 1:length(GC_all)){
    GC_prop_all[[i]] <- movingproportion(GC_all[[i]],movingwindow)
    GC_prop_all[[i]] <- GC_prop_all[[i]][seq(1,nrow(GC_prop_all[[i]]),10),]
    colnames(GC_prop_all[[i]]) <- c("chromosome","start","end","score")
    GC_prop_all[[i]] <- GRanges(seqnames = GC_prop_all[[i]]$chromosome,
                                ranges = IRanges(GC_prop_all[[i]]$start,
                                                 (GC_prop_all[[i]]$end + 9)),
                                score = GC_prop_all[[i]]$score)
    seqlengths(GC_prop_all[[i]]) <- end(getWholeChromosomes(GC_prop_all[[i]]))
    setTxtProgressBar(pb,i)
  }
  save(GC_prop_all,file = paste0("~/work/Julia/ecotypes_comparison/rdata/",ecotype,"_",movingwindow,"_all.rdata"))
} else{
  load(paste0("~/work/Julia/ecotypes_comparison/rdata/",ecotype,"_",movingwindow,"_all.rdata"))
}
# Computing moving iNOME signal -------------------------------------------
if(ecotype=="Col0"){
  control <- GC_prop_all[grep(pattern = "CE",x = names(GC_prop_all))][[1]]
  treatment <- GC_prop_all[grep(pattern = "ME",x = names(GC_prop_all))][[1]]
} else{
  control <- GC_prop_all[grep(pattern = "0h",x = names(GC_prop_all))][[1]]
  treatment <- GC_prop_all[grep(pattern = "12h",x = names(GC_prop_all))][[1]]
}

treatment$score <- treatment$score - control$score
treatment <- treatment[!is.na(treatment$score)]
treatment <- treatment[is.finite(treatment$score)]
# Exporting as bigwig -----------------------------------------------------
export(treatment,paste0("~/work/Julia/ecotypes_comparison/moving_average/",ecotype,"_","all","_",movingwindow,".bw"),"bw")
####-####
###################################################
#### Exporting methylation profiles as bigwigs ####
###################################################
####-####
# Settings ----------------------------------------------------------------
ecotype <- "Col0"
# Loading files -----------------------------------------------------------
if(ecotype=="Col0"){
  files <- dir(path = paste0("~/work/Julia/cx_root_vs_leaf/root/"),pattern = "12h")
  cx_report <- lapply(paste0("~/work/Julia/cx_root_vs_leaf/root/",files),function(x){
    x <- readBismark2(x)
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  })
  names(cx_report) <- gsub(pattern = "_CX_report.txt",replacement = "",files)
} else{
  files <- dir(path = paste0("~/work/ecotypes/NOME/CX_reports/",ecotype,"/"))
  cx_report <- lapply(paste0("~/work/ecotypes/NOME/CX_reports/",ecotype,"/",files),function(x){
    x <- readBismark2(x)
    x <- keepSeqlevels(x,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
  })
  names(cx_report) <- gsub(pattern = "_R_fp_bismark_bt2_pe.deduplicated.CX_report.txt",replacement = "",files)
}
# Loading Genome ----------------------------------------------------------
if(ecotype=="Col0"){
  genome <- import(paste0("~/work/ecotypes/annotations/",ecotype,"/genome/",dir(paste0(path = "~/work/ecotypes/annotations/",ecotype,"/genome/"),pattern = "fasta")))[1:5]
  names(genome) <- gsub(pattern = " CHROMOSOME dumped from ADB: Feb/3/09 16:9; last updated: 2009-02-02",
                        replacement = "",names(genome))
  names(genome) <- gsub(pattern = " CHROMOSOME dumped from ADB: Feb/3/09 16:10; last updated: 2009-02-02",
                        replacement = "",names(genome))
} else{
  genome <- import(paste0("~/work/ecotypes/annotations/",ecotype,"/genome/",dir(paste0(path = "~/work/ecotypes/annotations/",ecotype,"/genome/"),pattern = "fasta")))[1:5]
}

# Extracting GC -----------------------------------------------------------
GC <- lapply(cx_report,extractGC,genome = genome,contexts = "ALL")

# Computing methylation profiles ------------------------------------------
if(ecotype=="Col0"){
  control <- GC[grep(pattern = "CE",x = names(GC))][[1]]
  treatment <- GC[grep(pattern = "ME",x = names(GC))][[1]]
} else{
  control <- GC[grep(pattern = "0h",x = names(GC))][[1]]
  treatment <- GC[grep(pattern = "12h",x = names(GC))][[1]]
}

contexts <- names(treatment)
treatment_profiles <- vector("list",length(treatment))
control_profiles <- vector("list",length(control))
names(treatment_profiles) <- contexts
names(control_profiles) <- contexts

for(i in 1:length(treatment)){
  regions <- getWholeChromosomes(treatment[[i]])
  treatment_profiles[[i]] <- vector("list",length(regions))
  control_profiles[[i]] <- vector("list",length(regions))
  names(treatment_profiles[[i]]) <- seqlevels(regions)
  names(control_profiles[[i]]) <- seqlevels(regions)
  for(j in 1:length(regions)){
    treatment_profiles[[i]][[j]] <- computeMethylationProfile(treatment[[i]],region = regions[j],windowSize = 15,context = contexts[i])
    control_profiles[[i]][[j]] <- computeMethylationProfile(control[[i]],region = regions[j],windowSize = 15,context = contexts[i])
  }
}

for(i in 1:length(treatment_profiles)){
  for(j in 1:length(treatment_profiles[[i]])){
    treatment_profiles[[i]][[j]]$Proportion <- treatment_profiles[[i]][[j]]$Proportion - control_profiles[[i]][[j]]$Proportion
  }
}
# Exporting profiles as bigwigs -------------------------------------------
iNOME_CG <- treatment_profiles[[1]]
iNOME_CHG <- treatment_profiles[[2]]
iNOME_CHH <- treatment_profiles[[3]]

iNOME_CG <- c(iNOME_CG[[1]],iNOME_CG[[2]],iNOME_CG[[3]],iNOME_CG[[4]],iNOME_CG[[5]])
iNOME_CHG <- c(iNOME_CHG[[1]],iNOME_CHG[[2]],iNOME_CHG[[3]],iNOME_CHG[[4]],iNOME_CHG[[5]])
iNOME_CHH <- c(iNOME_CHH[[1]],iNOME_CHH[[2]],iNOME_CHH[[3]],iNOME_CHH[[4]],iNOME_CHH[[5]])

iNOME_CG <- iNOME_CG[!is.na(iNOME_CG$Proportion)]
iNOME_CHG <- iNOME_CHG[!is.na(iNOME_CHG$Proportion)]
iNOME_CHH <- iNOME_CHH[!is.na(iNOME_CHH$Proportion)]

iNOME_CG <- iNOME_CG[is.finite(iNOME_CG$Proportion)]
iNOME_CHG <- iNOME_CHG[is.finite(iNOME_CHG$Proportion)]
iNOME_CHH <- iNOME_CHH[is.finite(iNOME_CHH$Proportion)]

iNOME_CG <- GRanges(seqnames(iNOME_CG),IRanges(start(iNOME_CG),end(iNOME_CG)),score = iNOME_CG$Proportion)
iNOME_CHG <- GRanges(seqnames(iNOME_CHG),IRanges(start(iNOME_CHG),end(iNOME_CHG)),score = iNOME_CHG$Proportion)
iNOME_CHH <- GRanges(seqnames(iNOME_CHH),IRanges(start(iNOME_CHH),end(iNOME_CHH)),score = iNOME_CHH$Proportion)

seqlengths(iNOME_CG) <- end(getWholeChromosomes(iNOME_CG))
seqlengths(iNOME_CHG) <- end(getWholeChromosomes(iNOME_CHG))
seqlengths(iNOME_CHH) <- end(getWholeChromosomes(iNOME_CHH))

export(iNOME_CG,paste0("~/work/Julia/ecotypes_comparison/mnase/",ecotype,"_iNOME_CG.bw"),format = "bw")
export(iNOME_CHG,paste0("~/work/Julia/ecotypes_comparison/mnase/",ecotype,"_iNOME_CHG.bw"),format = "bw")
export(iNOME_CHH,paste0("~/work/Julia/ecotypes_comparison/mnase/",ecotype,"_iNOME_CHH.bw"),format = "bw")





####-####
######################################################################
#### Calculating gene distribution accross genome for each method ####
######################################################################
####-####
# Loading DMRs ------------------------------------------------------------
## Loading DMRs
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")
## Removing irrelevent data
DMRs <- root_12h_CHH_ME_CE$`0.3`$`700`
## Cleaning global environment
rm(root_12h_CHH_ME_CE)
DMRs <- DMRs %>%
  mutate(proportion = proportion2 - proportion1) %>%
  select(context,proportion)

# Loading ATAC peaks -------------------------------------------------------
## Loading ATAC-seq peaks
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak") %>%
  keepSeqlevels(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse") %>%
  select(score,signalValue) %>%
  unique()
# Loading MNase data ------------------------------------------------------
MNase_data <- import("~/work/Julia/ecotypes_comparison/mnase/mNASE_col_50.bw") %>%
  dropSeqlevels(.,c("Chrmitochondria","Chrchloroplast"),"coarse")

MNase_peaks <- import("/Volumes/PhD_one/PhD/data/mnase/data/SRR957780_q005_peaks.narrowPeak") %>%
  unique() %>%
  dropSeqlevels(c("chloroplast", "mitochondria"),"coarse") %>%
  convert_cx_data() %>%
  mutate(type = "MNase") %>%
  select(type,ID = name,!score:peak)
# Loading fasta -----------------------------------------------------------
## Loading fasta
Fasta <- import("~/work/ecotypes/annotations/Col0/TAIR10_chr_all.fasta")
names(Fasta) <- c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5")
Chromosome.lengths <- width(Fasta)
## Removing Chloroplast and Mitochondria data
Fasta <- Fasta[1:5]
# Loading GFF -------------------------------------------------------------
## Loading GFF
GFF <- import("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff") %>%
  filter(type == "gene") %>%
  select(type,ID) %>%
  dropSeqlevels(.,c("ChrC", "ChrM"),"coarse")
# Calculating average iNOME per gene --------------------------------------
iNOME_control <- readBismark("~/work/ecotypes/NOME/CX_reports/Col0/CE_12h_R.R1.RUN0076.L3._CX_report.txt") %>%
  dropSeqlevels(.,c("Mt", "Pt"),"coarse") %>%
  convert_cx_data(.,"CHH") %>%
  extract_GC(.,Fasta,Chromosome.lengths) %>%
  analyseReadsInsideRegionsForCondition(GFF,.,"CHH",cores = 4) %>%
  select(ID,proportionCHH)

iNOME_treatment <- readBismark("~/work/ecotypes/NOME/CX_reports/Col0/ME_12h_R.R1.RUN0076.L3._CX_report.txt") %>%
  dropSeqlevels(.,c("Mt", "Pt"),"coarse") %>%
  convert_cx_data(.,"CHH") %>%
  extract_GC(.,Fasta,Chromosome.lengths) %>%
  analyseReadsInsideRegionsForCondition(GFF,.,"CHH",cores = 4) %>%
  select(ID,proportionCHH)

average_gene_iNOME <- iNOME_treatment %>%
  mutate(iNOME = proportionCHH - iNOME_control$proportionCHH) %>%
  select(ID,iNOME) %>%
  as_tibble() %>%
  select(ID,iNOME)

# Calculating average MNase per gene --------------------------------------
overlaps <- findOverlaps(MNase_data,GFF)
mnase_data_list <- splitAsList(MNase_data[queryHits(overlaps)],subjectHits(overlaps))
pb <- txtProgressBar(0,length(mnase_data_list),style = 3)
average_gene_mnase <- GRanges()
for(i in 1:length(mnase_data_list)){
  average_gene_mnase <- c(average_gene_mnase,GRanges(seqnames(mnase_data_list[[i]][1]),
                                                     IRanges(min(start(mnase_data_list[[i]])),
                                                             max(end(mnase_data_list[[i]]))),
                                                     score = mean(mnase_data_list[[i]]$score)))
  setTxtProgressBar(pb,i)
}

average_gene_mnase$ID <- GFF$ID
average_gene_mnase <- average_gene_mnase %>%
  as_tibble() %>%
  select(ID,MNase = score)
# Classifying genes ------------------------------------------------------
## Overlapping DMRs with GFF
iNOME_genes <- GFF[overlapsAny(GFF,DMRs)]
ATAC_genes <- GFF[overlapsAny(GFF,ATAC_peaks)]

common_genes <- iNOME_genes[iNOME_genes$ID %in% ATAC_genes$ID] %>%
  mutate(type = "Common") %>%
  as_tibble() %>%
  inner_join(average_gene_iNOME,by = "ID") %>%
  inner_join(average_gene_mnase,by = "ID") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
unique_genes_iNOME <- iNOME_genes[!iNOME_genes$ID %in% ATAC_genes$ID] %>%
  mutate(type = "iNOME") %>%
  as_tibble() %>%
  inner_join(average_gene_iNOME,by = "ID") %>%
  inner_join(average_gene_mnase,by = "ID") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)
unique_genes_ATAC <- ATAC_genes[!ATAC_genes$ID %in% iNOME_genes$ID] %>%
  mutate(type = "ATAC") %>%
  as_tibble() %>%
  inner_join(average_gene_iNOME,by = "ID") %>%
  inner_join(average_gene_mnase,by = "ID") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

not_Identified <- GFF[!GFF$ID %in% iNOME_genes$ID]
not_Identified <- not_Identified[!not_Identified$ID %in% ATAC_genes$ID] %>%
  mutate(type = "Not Identified") %>%
  as_tibble() %>%
  inner_join(average_gene_iNOME,by = "ID") %>%
  inner_join(average_gene_mnase,by = "ID") %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

query <- c(common_genes,unique_genes_iNOME,unique_genes_ATAC,MNase_peaks,not_Identified)

query$type <- factor(query$type,c("Common","ATAC","iNOME","Not Identified","MNase"))

centromeres = tibble(chr = str_c("Chr",1:5), position = c(15086045, 3607929, 13587786, 3956021, 11725024))
heterochromatin = read_tsv("~/work/Javier/heterochomatin.bed", col_names = c("chr","start","end"), col_types = cols())
# Plotting gene density ---------------------------------------------------
query %>%
  as_tibble() %>%
  ggplot(aes(x = start/1e+6,colour = type,fill = type)) +
  facet_chr(chr = seqnames, rows = type) +
  geom_density(adjust = 1/10, alpha = 0.2) +
  scale_color_manual(values = c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4")) + 
  scale_fill_manual(values = c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4")) +
  geom_cytoband(chr_col = seqnames) +
  guides(color = FALSE, fill = FALSE)
ggsave("~/work/Figures/gene_density.pdf")

query %>%
  as_tibble() %>%
  filter(type != "MNase") %>%
  ggplot(aes(x = type,fill = type)) +
  geom_bar() +
  geom_text(stat='count', aes(label=..count..), vjust=-1) +
  scale_color_manual(values = c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4")) + 
  scale_fill_manual(values = c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4")) +
  xlab("") +
  ylab("") +
  theme(legend.position = "none")
ggsave("~/work/Figures/gene_count_bar.pdf")

# Plotting MNase around genes ---------------------------------------------
## Normalising around genes
all_genes <- c(common_genes,unique_genes_ATAC,unique_genes_iNOME)

MNase_gene_matrix <- list()
for(i in c("Common", "ATAC", "iNOME")){
  regions <- filter(all_genes,type == i) %>%
    as_tibble() %>%
    sample_n(2000) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  MNase_gene_matrix[[i]] <- normalizeToMatrix(MNase_data,regions,value_column = "score",
                                              mean_mode = "absolute",
                                              extend = 2000,w = 10,target_ratio = 0.4,
                                              background = NA,smooth = T)
}

## Creating line plots
MNase_gene_matrix_average <- list()
for(i in c("Common", "ATAC", "iNOME")){
  MNase_gene_matrix_average[[i]] <- tibble(coordinate = seq(10,6670,10),score = colMeans(MNase_gene_matrix[[i]],na.rm = T),type = i)
}
MNase_gene_matrix_average <- bind_rows(MNase_gene_matrix_average)
MNase_gene_matrix_average$type <- factor(MNase_gene_matrix_average$type,levels = c("Common","ATAC","iNOME"))
MNase_gene_matrix_average %>%
  ggplot(aes(x = coordinate,y = score)) +
  geom_line(aes(colour = type)) +
  scale_color_manual(values = c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4")) +
  geom_vline(xintercept = c(2000,4670),col = "black",lty = "dashed") +
  scale_x_continuous(breaks = c(0,2000,4670,6670),labels = c("-2kb","TSS","TES","+2kb"))
ggsave("~/work/Figures/MNase_around_classified_genes.pdf")

## Creating heatmaps
MNase_gene_heatmaps <- c()

colour <- c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4")
x <- 1
for(i in c("Common", "ATAC", "iNOME")){
  colour_function <- colorRamp2(c(0,0.075,0.15),c("white","#a9a9a9",colour[x]))
  MNase_gene_heatmaps <- MNase_gene_heatmaps + EnrichedHeatmap(MNase_gene_matrix[[i]],
                                                               name = i, na_col="black",col = colour_function,
                                                               column_title=i,
                                                               top_annotation = HeatmapAnnotation(enriched = anno_enriched()),
                                                               use_raster = TRUE,raster_device = "png")
  x <- x + 1
}
pdf("~/work/Figures/MNase_around_classified_genes_heatmap.pdf")
draw(MNase_gene_heatmaps)
dev.off()
# Plotting average accesssibility of genes --------------------------------
combinations <- combn(c("Common", "ATAC", "iNOME","Not Identified"),2)
comparisons <- list()
for(i in 1:ncol(combinations)){
  comparisons[[i]] <- dput(combinations[,i])
}

query %>%
  as_tibble() %>%
  filter(MNase < 1) 

x <- GFF %>%
  filter(ID == "AT3G41768")
y <- findOverlaps(MNase_data,x)

MNase_data[queryHits(y)]



query %>%
  as_tibble() %>%
  filter(type !="MNase") %>%
  filter(MNase < 1) %>%
  ggboxplot(x = "type",y = "MNase",fill = "type",ylab = "MNase signal",outliers.shape = NA) +
  scale_fill_manual(values = c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4")) +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75),labels = c(0,0.25,0.5,0.75)) +
  coord_cartesian(ylim = c(0,1.25)) +
  stat_compare_means(method = "anova",label.y = 1.25) +
  stat_compare_means(label = "p.signif",comparisons = comparisons,method = "t.test") +
  theme(axis.title.x = element_blank(),legend.position = "none")
ggsave("~/work/Figures/average_MNase_classified_genes.pdf")
  

# Plotting average MNase of genes -----------------------------------------
combinations <- combn(c("Common", "ATAC", "iNOME","Not Identified"),2)
comparisons <- list()
for(i in 1:ncol(combinations)){
  comparisons[[i]] <- dput(combinations[,i])
}

query %>%
  as_tibble() %>%
  filter(type !="MNase") %>%
  ggboxplot(x = "type",y = "iNOME",fill = "type",ylab = "iNOME signal") +
  scale_fill_manual(values = c("#e6194B","#3cb44b","#4363d8","#f58231","#911eb4")) +
  scale_y_continuous(breaks = c(0,0.25,0.5),labels = c(0,0.25,0.5)) +
  coord_cartesian(ylim = c(0,1)) +
  stat_compare_means(method = "anova",label.y = 1) +        # Add global annova p-value
  stat_compare_means(label = "p.signif",comparisons = comparisons,method = "t.test") +
  theme(axis.title.x = element_blank(),legend.position = "none")

####-####
#############################################
#### Extracting cluster regions from GFF ####
#############################################
####-####
# Loading GFFs ------------------------------------------------------------
setwd("~/work/Julia/ecotypes_comparison/GFFs/")
files <- dir(pattern = "edited.gff3")[2:5]
files <- c("C24_custom_edited.gtf",files)
gffs <- list()
for(i in files){
  gffs[[i]] <- import(i) %>%
  {if(i == "C24_custom_edited.gtf") mutate(.,ID = gene_id) else .} %>%
    select(type,ID) %>%
    {if(i != "C24_custom_edited.gtf") filter(.,type == "gene") else .} %>%
    keepSeqlevels(c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
}
names(gffs) <- c("C24", "Col0", "Cvi", "Ler", "WS2")

gffs[["C24"]] <- gffs[["C24"]] %>%
  as_tibble() %>%
  group_by(ID) %>%
  summarise(seqnames,start = min(start),end = max(end),strand,type = "gene") %>%
  ungroup() %>%
  makeGRangesFromDataFrame(keep.extra.columns = T) %>%
  unique() %>%
  select(type,ID)

# Storing clustered genes -------------------------------------------------
clusters <- list()
clusters[[1]] <- c("AT5G47950", "AT5G47980", "AT5G47990", "AT5G48000", "AT5G48010")
clusters[[2]] <- toupper(c("At5g42580", "At5g42590", "At5g42600"))
clusters[[3]] <- toupper(c("At4g15300", "At4g15310", "At4g15330", "At4g15340", 
                           "At4g15350", "At4g15360", "At4g15370", "At4g15380", "At4g15390", 
                           "At4g15400"))
names(clusters) <- paste0("Cluster_",seq(1,length(clusters),1))

# Extracting genes from gffs ----------------------------------------------
widths <- c(25000,17500,70000)
clustered_regions <- list()
for(i in 1:5){
  clustered_regions[[i]] <- list()
  for(j in 1:3){
    clustered_regions[[i]][[j]] <- gffs[[i]] %>%
      filter(ID %in% clusters[[j]]) %>%
      as_tibble() %>%
      group_by(seqnames) %>%
      summarise(start = min(start),end = max(end),midpoint = (start + end)/2,name = paste0("cluster_",j)) %>%
      mutate(start = midpoint - widths[j],end = midpoint + widths[j]) %>%
      makeGRangesFromDataFrame(keep.extra.columns = T) %>%
      select(name)
    clustered_regions[[i]][[j]]$name <- paste0(clustered_regions[[i]][[j]]$name,".",seq(1,length(clustered_regions[[i]][[j]]$name),1))
  }
  clustered_regions[[i]] <- do.call("c",clustered_regions[[i]])
}
names(clustered_regions) <- names(gffs)

# Exporting regions -------------------------------------------------------
for(i in 1:length(clustered_regions)){
  export(clustered_regions[[i]],con = paste0("~/work/Julia/ecotypes_comparison/GFFs/",names(clustered_regions[i]),"_extended_region.bed"),format = "bed")
}
 

combinations <- combn(c("Col0","Cvi","C24","WS2","Ler"),2)
comparisons <- list()
for(i in 1:ncol(combinations)){
  comparisons[[i]] <- dput(combinations[,i])
}

####-####
##################
#### Plotting ####
##################
####-####
# Size of features density plot (1) --------------------------------------------
## Extracting common DMRs
cDMRs <- intersect(DMRs,ATAC_peaks)
cDMRs <- intersect(cDMRs,DNase_peaks)
cDMRs$type <- "cDMR"
## Extracting unique DMRs
sDMRs <- setdiff(DMRs,ATAC_peaks)
sDMRs <- setdiff(sDMRs,DNase_peaks)
sDMRs$type <- "sDMR"
## Extracting common ATAC-seq peaks
cATAC <- intersect(ATAC_peaks,DMRs)
cATAC <- intersect(ATAC_peaks,DNase_peaks)
cATAC$type <- "cATAC"
## Combining data
feature_size <- as.data.frame(unlist(GRangesList(cDMRs,sDMRs,cATAC)))[c(4,6)]
## Plotting data
ggplot(feature_size,aes(x = width,colour = type)) +
  geom_density() +
  xlim(0,3000) +
  xlab("Size of features (bp)") +
  ylab("Density") +
  theme_cowplot()
ggsave("~/work/ecotypes/figures/density_plot.pdf")
# Annotated DMRs bar chart (1) ------------------------------------------------
## Unlisting GRangeList
cDMRs_anno <- unlist(cDMRs_anno)
sDMRs_anno <- unlist(sDMRs_anno)
## Counting number of DMRs per feature
DMRs_stats <- data.frame(Features = character(),Percentage = integer(),Type = character())
features <- unique(cDMRs_anno$type)
for(i in 1:length(unique(cDMRs_anno$type))){
  x <- data.frame(Features = features[i],Percentage = length(cDMRs_anno[cDMRs_anno$type==features[i]]),Type = "Common DMR")
  y <- data.frame(Features = features[i],Percentage = length(sDMRs_anno[sDMRs_anno$type==features[i]]),Type = "Specfic DMR")
  DMRs_stats <- rbind(DMRs_stats,x,y)
}
DMRs_stats[DMRs_stats$Type=="Common DMR",]$Percentage <- (DMRs_stats[DMRs_stats$Type=="Common DMR",]$Percentage / sum(DMRs_stats[DMRs_stats$Type=="Common DMR",]$Percentage)) * 100
DMRs_stats[DMRs_stats$Type=="Specfic DMR",]$Percentage <- (DMRs_stats[DMRs_stats$Type=="Specfic DMR",]$Percentage / sum(DMRs_stats[DMRs_stats$Type=="Specfic DMR",]$Percentage)) * 100
## Plotting data
ggplot(DMRs_stats,aes(x = Type, y = Percentage, fill = Features)) +
  geom_bar(stat = "identity",color = "black") +
  xlab("") +
  theme_cowplot()
ggsave("~/work/ecotypes/figures/DMR_comp_bar.pdf")
# Plotting metaplot (2) ----------------------------------------------------
## Converting matrix to dataframe
normalised_cDMRs_df <- data.frame()
for(i in 1:length(normalised_cDMRs)){
  x <- data.frame(score = (colMeans(normalised_cDMRs[[i]],na.rm = T)),
                  feature = unique(feature_scores[[i]]$type))
  x$position <- seq(10,(nrow(x)*10),10)
  x$type <- "Common DMR"
  normalised_cDMRs_df <- rbind(normalised_cDMRs_df,x)
}
normalised_sDMRs_df <- data.frame()
for(i in 1:length(normalised_sDMRs)){
  x <- data.frame(score = (colMeans(normalised_sDMRs[[i]],na.rm = T)),
                  feature = unique(feature_scores[[i]]$type))
  x$position <- seq(10,(nrow(x)*10),10)
  x$type <- "Specific DMR"
  normalised_sDMRs_df <- rbind(normalised_sDMRs_df,x)
}
normalised_df <- rbind(normalised_cDMRs_df,normalised_sDMRs_df)
x <- ggplot(normalised_df[normalised_df$feature=="DNase",],aes(x = position,y = score,colour = type)) +
  geom_line() +
  ggtitle("DNase") + 
  geom_vline(xintercept = 2000,linetype = "dashed") +
  scale_x_continuous(breaks = c(0,2000,4000),labels = c("-2kb","0bp","+2kb")) +
  theme_cowplot()
y <- ggplot(normalised_df[normalised_df$feature=="ATAC",],aes(x = position,y = score,colour = type)) +
  geom_line() +
  ggtitle("ATAC") + 
  geom_vline(xintercept = 2000,linetype = "dashed") +
  scale_x_continuous(breaks = c(0,2000,4000),labels = c("-2kb","0bp","+2kb")) +
  theme_cowplot()
f <- ggplot(normalised_df[normalised_df$feature=="H3K27me3",],aes(x = position,y = score,colour = type)) +
  geom_line() +
  ggtitle("H3K27me3") + 
  geom_vline(xintercept = 2000,linetype = "dashed") +
  scale_x_continuous(breaks = c(0,2000,4000),labels = c("-2kb","0bp","+2kb")) +
  theme_cowplot()
l <- ggplot(normalised_df[normalised_df$feature=="H3K27ac",],aes(x = position,y = score,colour = type)) +
  geom_line() +
  ggtitle("H3K27ac") + 
  geom_vline(xintercept = 2000,linetype = "dashed") +
  scale_x_continuous(breaks = c(0,2000,4000),labels = c("-2kb","0bp","+2kb")) +
  theme_cowplot()
ggarrange(x,y,f,l,nrow = 2,ncol = 2,common.legend = T,legend = "right")
ggsave("~/work/ecotypes/figures/MNase_metaplot.pdf")
# Plotting MCViP activity in plant lines (3) ----------------------------------
# Creating empty data frame
MCViP_df <- data.frame(Chromosome=factor(),
                       Proportion = numeric(),
                       Coordinate = numeric(),
                       Context = factor(),
                       Line = factor(),
                       Type = factor(),
                       stringsAsFactors=FALSE)
# Combining data frames to allow for easy plotting
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
# Removing Irrelevent data
MCViP_df <- MCViP_df[MCViP_df$Chromosome=="Chr1",]
MCViP_df <- MCViP_df[MCViP_df$Line!="LCol_0",]
# Plotting data
ggplot(MCViP_df,aes(x = Coordinate/1000000,y = Proportion,colour = Line)) +
  geom_line() +
  facet_grid(cols = vars(Context),rows = vars(Type)) +
  xlab("Genomic Coordinate (Mb)") +
  ylab("Methylation Proportion") +
  theme_cowplot()
ggsave("~/work/Figures/MCViP_activity_500_smooth.pdf",height = 10,width = 14)
# Plotting root vs leaf methylation profiles (4) ------------------------------
# Creating empty data frame to store data
root_CE_profile <- root_profile[grep(pattern = "CE",x = names(root_profile))]
root_ME_profile <- root_profile[grep(pattern = "ME",x = names(root_profile))]

leaf_CE_profile <- leaf_profile[grep(pattern = "CE",x = names(leaf_profile))]
leaf_ME_profile <- leaf_profile[grep(pattern = "ME",x = names(leaf_profile))]

root_iNOME_profile <- root_ME_profile
for(i in 1:length(root_iNOME_profile)){
  for(j in 1:length(root_iNOME_profile[[i]])){
    for(k in 1:length(root_iNOME_profile[[i]][[j]])){
      root_iNOME_profile[[i]][[j]][[k]]$Proportion <- root_iNOME_profile[[i]][[j]][[k]]$Proportion - root_CE_profile[[i]][[j]][[k]]$Proportion
    }
  }
}
names(root_iNOME_profile) <- gsub(pattern = "ME_",replacement = "iNOME_",x = names(root_iNOME_profile))

leaf_iNOME_profile <- leaf_ME_profile
for(i in 1:length(leaf_iNOME_profile)){
  for(j in 1:length(leaf_iNOME_profile[[i]])){
    for(k in 1:length(leaf_iNOME_profile[[i]][[j]])){
      leaf_iNOME_profile[[i]][[j]][[k]]$Proportion <- leaf_iNOME_profile[[i]][[j]][[k]]$Proportion - leaf_CE_profile[[i]][[j]][[k]]$Proportion
    }
  }
}
names(leaf_iNOME_profile) <- gsub(pattern = "ME_",replacement = "iNOME_",x = names(leaf_iNOME_profile))

root_vs_leaf_df <- data.frame(Chromosome=factor(),
                              Proportion = numeric(),
                              Coordinate = numeric(),
                              Context = factor(),
                              Timepoint = factor(),
                              Tissue = factor(),
                              stringsAsFactors=FALSE)
# Storing root data
for(i in 1:length(root_iNOME_profile)){
  for(j in 1:length(root_iNOME_profile[[i]])){
    for(k in 1:length(root_iNOME_profile[[i]][[j]])){
      root_vs_leaf_df <- rbind(root_vs_leaf_df,
                               data.frame(Chromosome = seqnames(root_iNOME_profile[[i]][[j]][[k]]),
                                          Coordinate = (start(root_iNOME_profile[[i]][[j]][[k]]) + end(root_iNOME_profile[[i]][[j]][[k]]))/2,
                                          Proportion = root_iNOME_profile[[i]][[j]][[k]]$Proportion,
                                          Context = root_iNOME_profile[[i]][[j]][[k]]$context,
                                          Timepoint = names(root_iNOME_profile)[i],
                                          Tissue = "Root"))
    }
  }
}
# Storing leaf data
for(i in 1:length(leaf_iNOME_profile)){
  for(j in 1:length(leaf_iNOME_profile[[i]])){
    for(k in 1:length(leaf_iNOME_profile[[i]][[j]])){
      root_vs_leaf_df <- rbind(root_vs_leaf_df,
                               data.frame(Chromosome = seqnames(leaf_iNOME_profile[[i]][[j]][[k]]),
                                          Coordinate = (start(leaf_iNOME_profile[[i]][[j]][[k]]) + end(leaf_iNOME_profile[[i]][[j]][[k]]))/2,
                                          Proportion = leaf_iNOME_profile[[i]][[j]][[k]]$Proportion,
                                          Context = leaf_iNOME_profile[[i]][[j]][[k]]$context,
                                          Timepoint = names(leaf_iNOME_profile)[i],
                                          Tissue = "Leaf"))
    }
  }
}

# Removing irrelevent data
root_vs_leaf_df <- root_vs_leaf_df[root_vs_leaf_df$Context=="CG",]
root_vs_leaf_df <- root_vs_leaf_df[root_vs_leaf_df$Chromosome=="Chr1",]
root_vs_leaf_df$Timepoint <- factor(root_vs_leaf_df$Timepoint,levels = c("iNOME_6h","iNOME_12h","iNOME_24h"))
# Cleaning global environment
rm(i,j,k)
# Plotting data
ggplot(root_vs_leaf_df,aes(x = Coordinate/1000000,y = Proportion,colour = Tissue)) +
  geom_line() +
  facet_grid(cols = vars(Timepoint),scales = "free") +
  xlab("Genomic Coordinate (Mb)") +
  ylab("Methylation Proportion") +
  theme_cowplot()
ggsave("~/work/Figures/root_vs_leaf_iNOME.pdf",height = 6,width = 14)
# Plotting Venn diagrams (5) --------------------------------------------------
pdf("~/work/Figures/Venn_root_CpG_pair.pdf",height = 10,width = 10)
draw.pairwise.venn(total_width_root_CpG[1,2],total_width_root_CpG[2,2],DMRs_ATAC_root_CpG,
                   category = c(paste0("iNOME","\n",total_width_root_CpG$label[1],"Mb"),
                                paste0("ATAC","\n",total_width_root_CpG$label[2],"Mb")),col = "transparent",
                   fill = c("orangered","steelblue"),alpha = 0.7)
dev.off()

pdf("~/work/Figures/Venn_root_CHG_pair.pdf",height = 10,width = 10)
draw.pairwise.venn(total_width_root_CHG[1,2],total_width_root_CHG[2,2],DMRs_ATAC_root_CHG,
                   category = c(paste0("iNOME","\n",total_width_root_CHG$label[1],"Mb"),
                                paste0("ATAC","\n",total_width_root_CHG$label[2],"Mb")),col = "transparent",
                   fill = c("orangered","steelblue"),alpha = 0.7)
dev.off()

pdf("~/work/Figures/Venn_root_CHH_pair.pdf",height = 10,width = 10)
draw.pairwise.venn(total_width_root_CHH[1,2],total_width_root_CHH[2,2],DMRs_ATAC_root_CHH,
                   category = c(paste0("iNOME","\n",total_width_root_CHH$label[1],"Mb"),
                                paste0("ATAC","\n",total_width_root_CHH$label[2],"Mb")),col = "transparent",
                   fill = c("orangered","steelblue"),alpha = 0.7)
dev.off()
# Plotting barplot for root and leaf annotation (6) -----------------------
## Combining root and leaf data frames
root_overlap$Tissue <- "Root"
leaf_overlap$Tissue <- "Leaf"
annotated_DMRs <- rbind(root_overlap[2:6,],leaf_overlap[2:6,])
ggplot(annotated_DMRs,aes(x = Tissue, y = Length/1000000, fill = Type)) +
  geom_bar(position="fill",stat="identity",colour = "black") +
  ylab("Accessible DNA (%)") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = seq(0,100,25)) +
  theme_cowplot()
ggsave("~/work/Figures/Annotated_root_leaf_DMRS.pdf",height = 8,width = 8)
# Plotting metaplot around genes (8) --------------------------------------
GC_df <- data.frame(Score = numeric(),Position = numeric(),
                    Sample = factor(),Context = factor())
for(i in 1:3){
  for(j in 1:3){
    GC_df <- rbind(GC_df,GC_mat[[i]][[j]])
  }
}
GC_df$Sample <- gsub(pattern = "ME_",replacement = "",x = GC_df$Sample)
GC_df$Sample <- factor(GC_df$Sample,levels = c("6h","12h","24h")) 

ggplot(GC_df,aes(x = Position, y = Score,colour = Context))+
  geom_line() +
  geom_vline(xintercept = c(GC_mat[[i]][[j]][rownames(GC_mat[[i]][[j]])=="d1",]$Position,
                            GC_mat[[i]][[j]][rownames(GC_mat[[i]][[j]])=="u100",]$Position),
             linetype = "dotted",colour = "red") +
  scale_x_continuous(breaks = c(0,1000,2340,3330),labels = c("-1kb","TSS","TES","+1kb")) +
  facet_grid(cols = vars(Sample)) +
  xlab("Genomic Coordinate") +
  ylab("iNOME-seq signal") +
  theme_cowplot()
ggsave("~/work/Figures/iNOME_genes_all.pdf",height = 6,width = 16)

# Plotting barplot for root and leaf annotation (9) -----------------------
## Combining root and leaf data frames
cDMRs_overlap$DMR <- "cDMRs"
sDMRs_overlap$DMR <- "sDMRs"
annotated_DMRs <- rbind(cDMRs_overlap[2:6,],sDMRs_overlap[2:6,])
ggplot(annotated_DMRs,aes(x = DMR, y = Length/1000000, fill = Type)) +
  geom_bar(position="fill",stat="identity",colour = "black") +
  ylab("Accessible DNA (%)") +
  scale_y_continuous(breaks = c(0,0.25,0.5,0.75,1),labels = seq(0,100,25)) +
  theme_cowplot()
ggsave("~/work/Figures/Annotated_cDMR_sDMR.pdf",height = 8,width = 8)
# Plotting boxplot of classified genes (10) -------------------------------
stat.test <- allClasses %>%
  tukey_hsd(tpm ~ class) %>%
  add_xy_position(x = "class") %>%
  mutate(y.position = seq(0.00009,0.00014,0.00001))
ggboxplot(allClasses,x = "class",y = "tpm",scales = "free",outlier.shape = NA) +
  scale_y_continuous(limits = c(0,0.00014)) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif",tip.length = 0.00005) +
  theme_cowplot() +
  xlab("Class") +
  ylab("TPM")
ggsave("~/iNOME/Figures/Main/Figure2/Parts/fig2G.pdf")
# Plotting cDMRs, sDMRs and all DMRs (11) -------------------------------------
ggplot(normalised_df,aes(x = Position,y = Score)) +
  geom_line() +
  facet_grid(cols = vars(Type),rows = vars(Feature),scales = "free") +
  geom_vline(xintercept = c(1000,1870),linetype = "dashed", colour = "red") +
  annotate("segment",x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size = 1.2) +
  scale_x_continuous(breaks = c(0,1435,2860),labels = c("-1kb","iNOME regions","+1kb")) +
  theme_cowplot() +
  theme(strip.background = element_blank(),axis.title.x = element_blank())
ggsave("~/work/Figures/features_around_DMRs.pdf",height = 10,width = 10)
# Plotting cDMRs, sDMRs and all DMRs (Intergenic)(12) ---------------------
ggplot(normalised_df,aes(x = Position,y = Score)) +
  geom_line() +
  facet_grid(cols = vars(Type),rows = vars(Feature),scales = "free") +
  geom_vline(xintercept = c(1000,1870),linetype = "dashed", colour = "red") +
  annotate("segment",x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size = 1.2) +
  scale_x_continuous(breaks = c(0,1435,2860),labels = c("-1kb","iNOME regions","+1kb")) +
  theme_cowplot() +
  theme(strip.background = element_blank(),axis.title.x = element_blank())
ggsave("~/work/Figures/features_around_DMRs_intergenic.pdf",height = 10,width = 10)


# Plotting iNOME around TSS heat map (13) ---------------------------------
iNOME_12h <- iNOME[grep("12h",names(iNOME))]
mat_12h <- vector("list",length(iNOME_12h))
iNOME_heatmap <- c()
pb <- txtProgressBar(0,length(iNOME_12h[[1]]),style = 3)
for(i in 1:length(iNOME_12h[[1]])){
  mat_12h[[i]] <- normalizeToMatrix(iNOME_12h[[1]][[i]],tss,value_column = "score",
                                    extend = c(500,1000), mean_mode = "absolute",
                                    w = 10)
  iNOME_heatmap <- iNOME_heatmap + EnrichedHeatmap(mat_12h[[i]],
                                                   col = colorRamp2(c(0,0.6), c("white","brown")),
                                                   name = paste0("G",names(iNOME_12h[[1]])[i]),
                                                   column_title = paste0("G",names(iNOME_12h[[1]])[i]),
                                                   row_order = order(as.numeric(tss$rank)))
  setTxtProgressBar(pb,i)
}

hm <- iNOME_heatmap + Heatmap(log2(tss$TPM+1), 
                              col = c("white", "orange"), 
                              name = "log2(TPM+1)", 
                              show_row_names = FALSE, width = unit(5, "mm"))

pdf("~/work/Figures/iNOME_TSS_100win.pdf",width = 10,height=10)
draw(hm,ht_gap = unit(1,"cm"))
dev.off()

# Plotting features around TSS (14) ---------------------------------------
ggplot(mat_12h_minus_df,aes(x = Position,y = Score + 0.5,colour = Type)) +
  geom_line() +
  geom_line(data = mat_mNASE_df,mapping = aes(x = Position,y = Score)) +
  scale_x_continuous(breaks = c(0,500,1500),labels = c("-0.5kb","TSS","+1kb")) +
  geom_vline(xintercept = 500,colour = "red",linetype = "dashed") +
  ylab("Score") + 
  theme_cowplot()
ggsave("~/work/Figures/iNOME_mNASE.pdf",height = 5,width = 5)
# Plotting barchart (15) --------------------------------------------------
ggplot(regulation_df,aes(x = Type, y = Value,fill = Regulation)) +
  geom_bar(position="fill", stat="identity") +
  theme_cowplot()
ggsave("~/work/Julia/Figures/Col_Ler_deg_barchart_0.001_700.pdf",height = 10,width = 10)
ggplot(fisher_pvalue,aes(x = Type,y = P.value)) +
  geom_point() +
  geom_hline(yintercept = 0.05,linetype = "dotted",col = "red") +
  theme_cowplot()
ggsave("~/work/Julia/Figures/Col_Ler_deg_fisher_0.001_700.pdf",height = 10,width = 10)
# Plotting Volcano plot (15) ----------------------------------------------
pdf("~/work/Col0_vs_Ler_volcano.pdf",height = 7,width = 14)
EnhancedVolcano(res,lab = rownames(res),x = 'log2FoldChange',y = 'pvalue',
                pCutoff = 0.001,title = "Col0 versus LER",subtitle = "",
                legendPosition = "right",pointSize = 0.5,colAlpha = 0.3)
dev.off()

EnhancedVolcano(res,lab = rownames(res),x = 'log2FoldChange',y = 'padj'
                pCutoff = 0.05,title = "WT versus RNAi",subtitle = "",
                legendPosition = "right",pointSize = 0.5,colAlpha = 0.3)
# Plotting comparison between iNOME and MNase (16) ------------------------
ggplot(iNOME_MNase_comparison,aes(x = Position,y = Score,colour = Type)) + 
  geom_line() + 
  stat_cor() +
  facet_grid(rows = vars(quartile)) +
  scale_x_continuous(breaks = c(0,501,1001,1501),labels = c("-0.5kb","TSS","+0.5kb","+1kb")) +
  scale_y_continuous("MNase",breaks = seq(0,0.2,0.05),labels = seq(0,0.2,0.05),sec.axis = sec_axis(~.+0.71,name = "1 - iNOME-seq",breaks = seq(0,0.9,0.1),labels = seq(0,0.9,0.1))) +
  scale_color_manual(values = c("firebrick","dodgerblue")) + 
  geom_vline(xintercept = 501,linetype = "dashed") +
  annotate("segment",x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size = 1.2) +
  theme_cowplot()
ggsave("~/work/Figures/iNOME_mnase_comp.pdf",height = 5,width  = 9)
####-####
