rm(list=ls())
source("~/iNOME/Scripts/R/libraries.R")
# Editing counts files ----------------------------------------------------
# Loading counts files
sampleTable <- read.table("/Volumes/shared/Shared210/Ryan/ecotypes_rna/RNA/counts/sampleTable.txt",header = T)
files <- paste0("/Volumes/shared/Shared210/Ryan/ecotypes_rna/RNA/counts/original/",sampleTable[,2])
ecotype_count <- lapply(files,read_delim,delim = "\t",col_names = c("gene","count"))
names(ecotype_count) <- paste0(sampleTable$condition,"_",seq(1,2,1))
# Loading blast tables
files <-  paste0("~/work/Julia/ecotypes_comparison/blast/",dir(path = "~/work/Julia/ecotypes_comparison/blast"))
blast_tables <- list()
for(i in 1:length(files)){
  x <- read_delim(files[i],delim = "\t",col_names = F) %>% 
    select(gene = X1,Col0 = X2,pident = X3)
  x <- x %>%
    group_split(Col0)
  for(j in 1:length(x)){
    x[[j]] <- x[[j]] %>%
      arrange(desc(pident)) %>%
      slice(n = 1)
  }
  blast_tables[[i]] <- bind_rows(x) %>%
    select(-pident)
}
names(blast_tables) <- c("Cvi","Ler","WS2")
blast_tables[["WS2"]] <- blast_tables[["WS2"]] %>%
  bind_rows(tibble(gene = "AT5G36130",Col0 = "AT5G36130"))
# Loading C24 genes
C24_genes <- import("~/work/Julia/ecotype_comparisons_hub/C24/Annotation/C24_ragtag.gff3") %>%
  as_tibble() %>%
  filter(type == "gene") %>%
  select(C24 = Name,Col0 = Name) %>%
  bind_rows(tibble(C24 = "AT5G47980",Col0 = "AT5G47980"))
# extracting common genes
blast_tables <- blast_tables[[1]] %>%
  inner_join(blast_tables[[2]],by = c("Col0" = "Col0")) %>%
  inner_join(blast_tables[[3]],by = c("Col0" = "Col0")) %>%
  inner_join(C24_genes,by = c("Col0" = "Col0")) %>%
  rename(Cvi = gene.x,Ler = gene.y,WS2 = gene) 
blast_tables <- list(C24 = tibble(gene = blast_tables$C24,Col0 = blast_tables$Col0),Cvi = tibble(gene = blast_tables$Cvi,Col0 = blast_tables$Col0),
                     Ler = tibble(gene = blast_tables$Ler,Col0 = blast_tables$Col0),WS2 = tibble(gene = blast_tables$WS2,Col0 = blast_tables$Col0),
                     Col0 = tibble(gene = blast_tables$Col0,Col0 = blast_tables$Col0))
# Converting gene names to Col0 gene names
files <- paste0("/Volumes/shared/Shared210/Ryan/ecotypes_rna/RNA/counts/edited/",sampleTable[,2])
for(i in 1:length(ecotype_count)){
  ecotype_blast <- blast_tables[grep(gsub("_.*","",names(ecotype_count)[i]),x = names(blast_tables))]
  x <- ecotype_count[[i]] %>%
    inner_join(ecotype_blast[[1]]) %>%
    mutate(gene = Col0) %>%
    select(gene,count)
  if(grepl("C24",names(ecotype_count)[i])){
    x <- bind_rows(x,tibble(gene = "AT5G47980",count = 0))
  }
  if(grepl("WS2",names(ecotype_count)[i])){
    x <- bind_rows(x,tibble(gene = "AT5G36130",count = 0))
  }
  x %>%
    arrange(gene) %>%
  write.table(file = files[i],append = F,quote = F,row.names = F,col.names = F,sep = "\t")
}

# DEseq2 ------------------------------------------------------------------
rm(list=ls())
# Loading data
sampleTable <- read.table("/Volumes/shared/Shared210/Ryan/ecotypes_rna/RNA/counts/sampleTable.txt",header = T)
sampleTable$sampleName <- paste0(sampleTable$condition,"_",seq(1,2,1))
#rownames(sampleTable) <- paste0(sampleTable$condition,"_",sampleTable$sampleName)

ddsHTSeq<-DESeqDataSetFromHTSeqCount(sampleTable=sampleTable,
                                     directory = "/Volumes/shared/Shared210/Ryan/ecotypes_rna/RNA/counts/edited/",
                                     design = ~ condition)
dds<-DESeq(ddsHTSeq)
# Comparing replicates to eachother
countsnorm<-counts(dds,normalized=TRUE)
## XY plots
Col0 <- ggplot(as_tibble(countsnorm),aes(x = Col0_1,y = Col0_2)) +
  geom_point() +
  geom_abline(intercept =0 , slope = 1,colour = "red",linetype = "dashed") +
  theme_cowplot()
C24 <- ggplot(as_tibble(countsnorm),aes(x = C24_1,y = C24_2)) +
  geom_point() +
  geom_abline(intercept =0 , slope = 1,colour = "red",linetype = "dashed") +
  theme_cowplot()
Cvi <- ggplot(as_tibble(countsnorm),aes(x = Cvi_1,y = Cvi_2)) +
  geom_point() +
  geom_abline(intercept =0 , slope = 1,colour = "red",linetype = "dashed") +
  theme_cowplot()
Ler <- ggplot(as_tibble(countsnorm),aes(x = Ler_1,y = Ler_2)) +
  geom_point() +
  geom_abline(intercept =0 , slope = 1,colour = "red",linetype = "dashed") +
  theme_cowplot()
WS2 <- ggplot(as_tibble(countsnorm),aes(x = WS2_1,y = WS2_2)) +
  geom_point() +
  geom_abline(intercept =0 , slope = 1,colour = "red",linetype = "dashed") +
  theme_cowplot()

ggarrange(Col0,C24,Cvi,Ler,WS2,ncol = 3,nrow = 2)
ggsave("~/iNOME/Figures/ecotypes_RNA/replicate_XY.pdf")
## PCA
rld <-rlog(dds)
pca <-prcomp(t(assay(rld)))
pov <- summary(pca)
pov <- pov$importance[2,]*100

PC12 <- as_tibble(pca$x) %>%
  mutate(sample = c("C24","C24","Cvi","Cvi","Ler","Ler","WS2","WS2","Col0","Col0")) %>%
  ggplot(aes(x = PC1,y = PC2,colour = sample)) +
  geom_point() +
  xlab(paste0("PC1 (",pov[1],"%)")) +
  ylab(paste0("PC2 (",pov[2],"%)")) +
  theme_cowplot()
PC13 <- as_tibble(pca$x) %>%
  mutate(sample = c("C24","C24","Cvi","Cvi","Ler","Ler","WS2","WS2","Col0","Col0")) %>%
  ggplot(aes(x = PC1,y = PC3,colour = sample)) +
  geom_point() +
  xlab(paste0("PC1 (",pov[1],"%)")) +
  ylab(paste0("PC3 (",pov[3],"%)")) +
  theme_cowplot()
PC23 <- as_tibble(pca$x) %>%
  mutate(sample = c("C24","C24","Cvi","Cvi","Ler","Ler","WS2","WS2","Col0","Col0")) %>%
  ggplot(aes(x = PC2,y = PC3,colour = sample)) +
  geom_point() +
  xlab(paste0("PC2 (",pov[2],"%)")) +
  ylab(paste0("PC3 (",pov[3],"%)")) +
  theme_cowplot()
ggarrange(PC12,PC13,PC23,common.legend = T,legend = "right")
ggsave("~/iNOME/Figures/ecotypes_RNA/replicate_PCA.pdf",height = 5,width = 10)
# Computing degs
## Creating vector of cluster genes
clusters <- list()
clusters[[1]] <- c("AT5G47950", "AT5G47980", "AT5G47990", "AT5G48000", "AT5G48010")
clusters[[2]] <- toupper(c("At5g42580", "At5g42590", "At5g42600"))
clusters[[3]] <- toupper(c("At5g36110", "At5g36130", "At5g36140", "At5g36150"))
clusters[[4]] <- toupper(c("At4g15300", "At4g15310", "At4g15330", "At4g15340", 
                           "At4g15350", "At4g15360", "At4g15370", "At4g15380", "At4g15390", 
                           "At4g15400"))
keyGenes <- unlist(clusters)
## Creating colours based on type
reference = "Col0" # "Col0" "C24" "Cvi" "Ler" "WS2"
treatment = "C24" # "Col0" "C24" "Cvi" "Ler" "WS2"
res = DESeq2::results(dds, contrast=c("condition",treatment,reference))
keyvals <- ifelse(rownames(res)%in%keyGenes,"red","grey")
names(keyvals)[keyvals == "red"] <- "cluster"
names(keyvals)[keyvals == "grey"] <- "non_cluster"
## Plotting Volcano plot

Volc1 <- EnhancedVolcano(res,lab = rownames(res),selectLab = c(keyGenes),x = "log2FoldChange",y = "pvalue",
                legendPosition = "right",title = c(paste0(reference," versus ",treatment)),labFace = "bold",
                pointSize = c(ifelse(rownames(res)%in%keyGenes,4,1)),colCustom = keyvals)

res = results(dds, contrast=c("condition","Cvi","Col0"))


ggsave(Volc1,filename = paste0("~/iNOME/Figures/ecotypes_RNA/",reference,"_",treatment,"_volcano.pdf"))

## Bar plots
reference = "WS2" # "Col0" "C24" "Cvi" "Ler" "WS2"
treatment = "Ler" # "Col0" "C24" "Cvi" "Ler" "WS2"
res = DESeq2::results(dds, contrast=c("condition",treatment,reference))
res_cluster_genes <- res[rownames(res) %in% keyGenes,]
res_cluster_genes <- tibble(genes = rownames(res_cluster_genes),fold_change = res_cluster_genes$log2FoldChange,padj = res_cluster_genes$padj,reference = reference,treatment = treatment) %>%
  mutate(fold_change = replace_na(fold_change,0),fold_change = ifelse(genes == "AT5G47980",0,fold_change),padj = ifelse(genes == "AT5G47980",1,padj),padj = replace_na(padj,1),Padj = ifelse(padj<0.05,"Significant","Not Significant"))
res_cluster_genes$cluster <- ifelse(res_cluster_genes$genes%in%clusters[[1]],"Thalianol",
                                    ifelse(res_cluster_genes$genes%in%clusters[[2]],"Marneral",
                                           ifelse(res_cluster_genes$genes%in%clusters[[3]],"Tirucalladienol",
                                                  ifelse(res_cluster_genes$genes%in%clusters[[4]],"Arabidiol","none"))))
res_cluster_genes$cluster <- factor(res_cluster_genes$cluster,levels = c("Thalianol","Marneral","Tirucalladienol","Arabidiol"))

ggplot(res_cluster_genes,aes(x = genes,y = fold_change,fill = Padj)) +
  ggtitle(paste0(reference," vs ",treatment)) +
  geom_bar(stat = "identity",size = 0.1) +
  annotate("segment",x = -Inf,xend = Inf,y = 0,yend = 0,size = 0.1) +
  scale_fill_manual(values = c("black","darkred")) +
  ylab("log2(fold change)") +
  xlab("Genes") +
  facet_wrap(vars(cluster),scales = "free_x") +
  theme_cowplot() +
  theme(strip.background = element_blank(),axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5))
ggsave(paste0("~/iNOME/Figures/ecotypes_RNA/","log_foldchange_",reference,"_",treatment,".pdf"))

comparisons <- tibble(reference = sort(rep(c("Col0","C24","Cvi","Ler","WS2"),5)),treatment = rep(c("Col0","C24","Cvi","Ler","WS2"),5)) %>%
  filter(reference != treatment)
all_clusters_ecotypes <- tibble()
for(i in 1:nrow(comparisons)){
  reference = comparisons %>%
    slice(i) %>%
    pull(reference)
  treatment = comparisons %>%
    slice(i) %>%
    pull(treatment)
  res = DESeq2::results(dds, contrast=c("condition",treatment,reference))
  res_cluster_genes <- res[rownames(res) %in% keyGenes,]
  res_cluster_genes <- tibble(genes = rownames(res_cluster_genes),fold_change = res_cluster_genes$log2FoldChange,padj = res_cluster_genes$padj,reference = reference,treatment = treatment) %>%
    mutate(fold_change = replace_na(fold_change,0),fold_change = ifelse(genes == "AT5G47980",0,fold_change),padj = ifelse(genes == "AT5G47980",1,padj),padj = replace_na(padj,1),Padj = ifelse(padj<0.05,"Significant","Not Significant"))
  res_cluster_genes$cluster <- ifelse(res_cluster_genes$genes%in%clusters[[1]],"Thalianol",
                                      ifelse(res_cluster_genes$genes%in%clusters[[2]],"Marneral",
                                             ifelse(res_cluster_genes$genes%in%clusters[[3]],"Tirucalladienol",
                                                    ifelse(res_cluster_genes$genes%in%clusters[[4]],"Arabidiol","none"))))
  res_cluster_genes$cluster <- factor(res_cluster_genes$cluster,levels = c("Thalianol","Marneral","Tirucalladienol","Arabidiol"))
  all_clusters_ecotypes <- bind_rows(all_clusters_ecotypes,res_cluster_genes)
}
all_clusters_ecotypes$cluster <- factor(all_clusters_ecotypes$cluster,levels = c("Thalianol","Marneral","Tirucalladienol","Arabidiol"))

for(i in c("Col0","C24","Cvi","Ler","WS2")){
  all_clusters_ecotypes %>%
    filter(reference == i) %>%
    ggplot(aes(x = genes,y = fold_change,fill = Padj)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values = c("black","darkred")) +
    ylab("log2(fold change)") +
    xlab("Genes") +
    facet_grid(cols = vars(cluster),rows = vars(treatment),scales = "free_x") + 
    geom_rect(aes(xmin = -Inf,xmax = Inf,ymin = -Inf,ymax = Inf),fill = NA,colour = "black",size = 0.1) +
    annotate("segment",x = -Inf,xend = Inf,y = 0,yend = 0,size = 0.1) +
    theme_cowplot() +
    theme(strip.background = element_blank(),strip.text.y = element_text(angle = 0),axis.text.x = element_text(angle = 45,vjust = 0.5,hjust = 0.5),
          axis.line = element_line(size = 0.1),axis.ticks = element_line(size = 0.1))
  ggsave(paste0("~/iNOME/Figures/ecotypes_RNA/log_foldchange_",i,"_all.pdf"),width = 14)
}
####-####

