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

.replace_chrs <- function(x){
  seqlevels(x) <- paste0("Chr",seqlevels(x))
  return(x)
}

unique_overlap_ranges <- function(gff,regions){
  average_DMR_proportion <- find_overlaps(gff,regions) %>%
    group_by(ID) %>%
    summarise(score = mean(proportion))
  gff <- inner_join(as_tibble(gff),as_tibble(average_DMR_proportion),by = "ID") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  return(gff)
}

extract_seq <- function(methylation.data,fasta,lengths,result = "GC"){
  output <- as_tibble(methylation.data) %>%
    mutate(start = ifelse(strand == "+",start - 1,start),end = ifelse(strand == "-",end + 1,end)) %>%
    group_split(seqnames) %>%
    .remove_bad_cytosines(.,lengths) %>%
    bind_rows() %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    .Extract_sequence(.,fasta) %>%
    filter(if(result == "GC") sequence == "GC" else TRUE) %>%
    mutate(start = ifelse(strand == "+",start + 1,start),end = ifelse(strand == "-",end - 1,end))
  return(output)
}

as_tibble(x) %>%
  mutate(start = ifelse(strand == "+",start - 1,start),end = ifelse(strand == "-",end + 1,end)) %>%
  group_split(seqnames) %>%
  .remove_bad_cytosines(.,lengths)
