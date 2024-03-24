rm(list=ls())
source("~/iNOME/Scripts/R/libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
# Loading Annotations -----------------------------------------------------
files <- dir("~/work/ecotypes/annotations/blast_edit",full.names = T)
gffs <- list()
for(i in files){
  gffs[[i]] <- import(i) %>%
    keepSeqlevels(.,paste0("Chr",seq(1,5,1)),"coarse")
}
files <- gsub("/Users/u1896702/work/ecotypes/annotations/blast_edit/","",files)
files <- gsub("*_.*","",files)
files[4] <- "Col0"
names(gffs) <- files

gffs_access <- gffs
for(i in 1:length(gffs_access)){
  end(gffs_access[[i]][strand(gffs_access[[i]])=="+"]) <- start(gffs_access[[i]][strand(gffs_access[[i]])=="+"]) + 150
  start(gffs_access[[i]][strand(gffs_access[[i]])=="+"]) <- start(gffs_access[[i]][strand(gffs_access[[i]])=="+"]) - 150
  start(gffs_access[[i]][strand(gffs_access[[i]])=="-"]) <- end(gffs_access[[i]][strand(gffs_access[[i]])=="-"]) - 150
  end(gffs_access[[i]][strand(gffs_access[[i]])=="-"]) <- end(gffs_access[[i]][strand(gffs_access[[i]])=="-"]) + 150
}
# Loading fasta files
files <- dir("~/work/ecotypes/annotations/fasta",full.names = T)
files <- files[-grep("C24",files)]
genomes <- lapply(files,function(x){ 
  x <- import(x)[1:5]
  names(x) <- paste0("Chr",seq(1,5,1))
  x
})
files <- gsub(".*/","",files)
names(genomes) <- gsub("\\..*","",files)
chromosome_lengths <- lapply(genomes,width)

# Extracting GC -----------------------------------------------------------
files <- dir("/Volumes/shared/Shared210/Ryan/ecotypes_nome/CX_reports",full.names = T)
control <- files[grep("0h",files)]
treatment <- c(files[grep("Col0_6h",files)],files[grep("12h",files)])
files <- c(control,treatment)
col0_files <- files[grep("Col0",files)]
files <- files[grep("Col0",files,invert = T)]
rm(control,treatment)
for(i in 1:length(col0_files)){
  Col0_data <- read_delim(col0_files[[i]],delim = "\t",col_names = c("seqnames","start","strand","readsM","readsN","context","Trinucleotide_context")) %>%
    filter(!is.na(seqnames)) %>%
    mutate(seqnames = paste0("Chr",seqnames)) %>%
    mutate(end = start) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    keepSeqlevels(.,value = paste0("Chr",seq(1,5,1)),pruning.mode = "coarse") %>%
    extract_seq(.,genomes[["Col0"]],chromosome_lengths[["Col0"]],"GC") %>%
    select(-sequence) %>%
    as_tibble()
  x <- gsub("_R_fp_bismark_bt2_pe.deduplicated.CX_report.txt","",col0_files)
  x <- gsub("/Volumes/shared/Shared210/Ryan/ecotypes_nome/CX_reports/MCViP_","",x)
  write_delim(Col0_data,paste0("/Volumes/shared/Shared210/Ryan/ecotypes_nome/GC/",x[i],"_GC.txt"),delim = "\t",col_names = F,append = F,quote_escape = F)
}
rm(Col0_data,col0_files)
for(i in 1:length(files)){
  genome <- str_extract(pattern = names(genomes),string = files[i])
  genome <- genome[!is.na(genome)]
  GC_data <- read_delim(files[[i]],delim = "\t",col_names = c("seqnames","start","strand","readsM","readsN","context","Trinucleotide_context")) %>%
    mutate(end = start) %>%
    makeGRangesFromDataFrame(keep.extra.columns = T) %>%
    keepSeqlevels(.,value = paste0("Chr",seq(1,5,1)),pruning.mode = "coarse") %>%
    extract_seq(.,genomes[[genome]],chromosome_lengths[[genome]],"GC") %>%
    select(-sequence) %>%
    as_tibble()
  x <- gsub("_R_fp_bismark_bt2_pe.deduplicated.CX_report.txt","",files)
  x <- gsub("/Volumes/shared/Shared210/Ryan/ecotypes_nome/CX_reports/MCViP_","",x)
  write_delim(GC_data,paste0("/Volumes/shared/Shared210/Ryan/ecotypes_nome/GC/",x[i],"_GC.txt"),delim = "\t",col_names = F,append = F,quote_escape = F)
}
rm(GC_data,files)

# Calculating avverage access ---------------------------------------------
files <- dir("/Volumes/shared/Shared210/Ryan/ecotypes_nome/GC",full.names = T)
GC_data <- list()
for(i in 1:length(files)){
  ecotype <- gsub("/Volumes/shared/Shared210/Ryan/ecotypes_nome/GC/","",files[i])
  ecotype <- gsub("_.*","",ecotype)
  x <- read_delim(files[i],delim = "\t",col_names = c("seqnames","start","end","width",
                                                      "strand","readsM","readsN","context",
                                                      "trinucleotide_context","sequence")) %>%
    select(-width,-sequence) %>%
    mutate(context = "CHH") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  genes <- gffs[[ecotype]] %>%
    filter(type == "gene")
  GC_data[[i]] <- analyseReadsInsideRegionsForCondition(genes,x,"CHH",cores = 4) %>%
    select(Name,score = proportionCHH)
}
files <- gsub("/Volumes/shared/Shared210/Ryan/ecotypes_nome/GC/","",files)
names(GC_data) <- gsub("_GC.txt","",files)
save(GC_data,file = "iNOME/Data/RData/average_meth_genes.RData")

control <- GC_data[grep("0h",names(GC_data))]
treatment <- GC_data[grep("0h",names(GC_data),invert = T)]
for(i in 1:length(treatment)){
  treatment[[i]]$score <- treatment[[i]]$score - control[[i]]$score
}
names(treatment) <- gsub("_.*","",names(treatment))
save(treatment,file = "iNOME/Data/RData/average_access_genes.RData")

files <- dir("/Volumes/shared/Shared210/Ryan/ecotypes_nome/GC",full.names = T)
GC_data <- list()
for(i in 1:length(files)){
  ecotype <- gsub("/Volumes/shared/Shared210/Ryan/ecotypes_nome/GC/","",files[i])
  ecotype <- gsub("_.*","",ecotype)
  x <- read_delim(files[i],delim = "\t",col_names = c("seqnames","start","end","width",
                                                      "strand","readsM","readsN","context",
                                                      "trinucleotide_context","sequence")) %>%
    select(-width,-sequence) %>%
    mutate(context = "CHH") %>%
    makeGRangesFromDataFrame(keep.extra.columns = T)
  genes <- gffs_access[[ecotype]] %>%
    filter(type == "gene")
  GC_data[[i]] <- analyseReadsInsideRegionsForCondition(genes,x,"CHH",cores = 4) %>%
    select(Name,score = proportionCHH)
}
files <- gsub("/Volumes/shared/Shared210/Ryan/ecotypes_nome/GC/","",files)
names(GC_data) <- gsub("_GC.txt","",files)
save(GC_data,file = "iNOME/Data/RData/average_meth_TSS.RData")

control <- GC_data[grep("0h",names(GC_data))]
treatment <- GC_data[grep("0h",names(GC_data),invert = T)]
for(i in 1:length(treatment)){
  treatment[[i]]$score <- treatment[[i]]$score - control[[i]]$score
}
names(treatment) <- gsub("_.*","",names(treatment))
save(treatment,file = "iNOME/Data/RData/average_access_tss.RData")

# Loading average accessibility -------------------------------------------
ecotypes <- c("Col0","Ler","Cvi","C24","WS2")
load("iNOME/Data/RData/average_access_genes.RData")
gene_access <- treatment
for(i in 1:length(gene_access)){
  gene_access[[i]] <- gene_access[[i]] %>%
    as_tibble() %>%
    select(id = Name,gene_access = score)
}
gene_access <- gene_access[ecotypes]

load("iNOME/Data/RData/average_access_tss.RData")
tss_access <- treatment
for(i in 1:length(tss_access)){
  tss_access[[i]] <- tss_access[[i]] %>%
    as_tibble() %>%
    select(id = Name,tss_access = score)
}
tss_access <- tss_access[ecotypes]
# Loading RNA data --------------------------------------------------------
ecotypes <- c("Col0","Ler","Cvi","C24","WS2")
RNA <- list()
for(i in ecotypes){
  file <- dir(paste0("/Volumes/shared/Shared210/Ryan/ecotypes_rna/RNA/counts/edited"),pattern = paste0(i,"_1_"),full.names = T)
  RNA[[i]] <- read_delim(file,delim = "\t",col_names = c("id","count"))
}
# Calculating TPM ---------------------------------------------------------
for(i in c(2,3,5)){
  gffs[[i]] <- gffs[[i]] %>%
    mutate(Parent = Name)
}

lengths <- list()
for(i in 1:length(gffs)){
  x <- gffs[[i]] %>%
    filter(type == "exon") %>%
    mutate(x = gsub("\\..*","",Parent),id = if(i %in% c(1,4)) unlist(x) else x) %>%
    select(-x)
  lengths[[i]] <- x %>%
    group_by(id) %>%
    group_split()
  lengths[[i]] <- lapply(lengths[[i]],function(x){
    x <- tibble(id = x$id[1],length = sum(width(reduce_ranges(x))))
    x
  })
  lengths[[i]] <- bind_rows(lengths[[i]])
}
names(lengths) <- names(gffs)
lengths <- lengths[ecotypes]

# Combining accessibility with expression ---------------------------------
for(i in 1:length(RNA)){
  RNA[[i]] <- RNA[[i]] %>%
    inner_join(lengths[[i]]) %>%
    inner_join(gene_access[[i]]) %>%
    inner_join(tss_access[[i]]) %>%
    mutate(tpm = count/length,tpm = tpm/(sum(tpm)/1000000)) %>%
    select(-count,-length) %>%
    filter_at(vars(gene_access,tss_access,tpm),all_vars(is.finite(.))) %>%
    mutate(gene_quartile = ntile(gene_access,3),tss_quartile = ntile(tss_access,3)) %>%
    mutate(gene_quartile = ifelse(gene_quartile == 3, "High",ifelse(gene_quartile == 2,"Medium",ifelse(gene_quartile == 1,"Low",gene_quartile))),
           tss_quartile = ifelse(tss_quartile == 3, "High",ifelse(tss_quartile == 2,"Medium",ifelse(tss_quartile == 1,"Low",tss_quartile))),
           ecotype = names(RNA)[i]) %>%
    filter(tpm < mean(tpm))
}
expression_accessibility <- bind_rows(RNA)
expression_accessibility$gene_quartile <- factor(expression_accessibility$gene_quartile,levels = c("Low","Medium","High"))
expression_accessibility$tss_quartile <- factor(expression_accessibility$tss_quartile,levels = c("Low","Medium","High"))
expression_accessibility$ecotype <- factor(expression_accessibility$ecotype,levels = ecotypes)
# Plotting comparison -----------------------------------------------------
stat.test <- expression_accessibility %>%
  group_by(ecotype) %>%
  tukey_hsd(tpm ~ tss_quartile) %>%
  add_xy_position(x = "ecotype") %>%
  mutate(y.position = rep(seq(185,201,7),5))

ggboxplot(expression_accessibility, x = "ecotype", y = "tpm", fill = "tss_quartile",
          scales = "free",outlier.shape = NA,size = 0.1) +
  scale_y_continuous(breaks = seq(0,180,20),labels = seq(0,180,20),limits = c(0,210)) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif",tip.length = 0.005,bracket.size = 0.1) +
  scale_fill_aaas() +
  labs(fill = "Accessibility") +
  xlab("Ecotype") +
  ylab("Expression") +
  theme_cowplot() +
  theme(axis.line = element_line(size = 0.1),axis.ticks = element_line(size = 0.1))
ggsave("iNOME/Figures/Main/Figure4/Parts/tss_expression_accessibility_comp.pdf")

stat.test <- expression_accessibility %>%
  group_by(ecotype) %>%
  tukey_hsd(tpm ~ gene_quartile) %>%
  add_xy_position(x = "ecotype") %>%
  mutate(y.position = rep(seq(185,201,7),5))

ggboxplot(expression_accessibility, x = "ecotype", y = "tpm", fill = "gene_quartile",
          scales = "free",outlier.shape = NA) +
  scale_y_continuous(breaks = seq(0,180,20),labels = seq(0,180,20),limits = c(0,210)) +
  stat_pvalue_manual(stat.test,label = "p.adj.signif",tip.length = 0.005,bracket.size = 0.1) +
  scale_fill_aaas() +
  labs(fill = "Accessibility") +
  xlab("Ecotype") +
  ylab("Expression") +
  theme_cowplot() +
  theme(axis.line = element_line(size = 0.1),axis.ticks = element_line(size = 0.1))
ggsave("iNOME/Figures/Main/Figure4/Parts/gene_expression_accessibility_comp.pdf")
