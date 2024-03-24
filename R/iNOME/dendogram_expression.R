rm(list=ls())
source("~/iNOME/Scripts/R/Libraries.R")
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
# Extracting keygenes -----------------------------------------------------
clusters <- list()
clusters[[1]] <- c("AT5G47950", "AT5G47980", "AT5G47990", "AT5G48000", "AT5G48010")
clusters[[2]] <- toupper(c("At5g42580", "At5g42590", "At5g42600"))
clusters[[3]] <- toupper(c("At5g36110", "At5g36130", "At5g36140", "At5g36150"))
clusters[[4]] <- toupper(c("At4g15300", "At4g15310", "At4g15330", "At4g15340", 
                           "At4g15350", "At4g15360", "At4g15370", "At4g15380", "At4g15390", 
                           "At4g15400"))
keyGenes <- unlist(clusters)
# Loading RNA data --------------------------------------------------------
ecotypes <- c("Col0","Ler","Cvi","C24","WS2")
RNA <- list()
for(i in ecotypes){
  file <- dir(paste0("/Volumes/Shared377/Ryan//iNOME/Ecotypes/ecotypes_rna/RNA/counts/edited"),pattern = paste0(i,"_1_"),full.names = T)
  RNA[[i]] <- read_delim(file,delim = "\t",col_names = c("id","count"))
}
for(i in c(2,3,5)){
  gffs[[i]] <- gffs[[i]] %>%
    mutate(Parent = Name)
}
# Extracting gene lengths -------------------------------------------------
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
# Calculating TPM ---------------------------------------------------------
for(i in 1:length(lengths)){
  RNA[[i]] <- RNA[[i]] %>%
    inner_join(lengths[[i]]) %>%
    mutate(tpm = count/length,tpm = tpm/(sum(tpm)/1000000)) %>%
    select(id,tpm)
}
all_expression <- RNA[[1]] %>%
  full_join(RNA[[2]],by = c("id" = "id")) %>%
  full_join(RNA[[3]],by = c("id" = "id")) %>%
  full_join(RNA[[4]],by = c("id" = "id")) %>%
  full_join(RNA[[5]],by = c("id" = "id"))
colnames(all_expression) <- c("id","Col0", "Ler", "Cvi", "C24", "WS2")


# Extracting keyGenes -----------------------------------------------------
x <- all_expression %>%
  filter(id %in% keyGenes) %>%
  as.data.frame(x)
x <- x %>%
  mutate(cluster = NA,cluster = ifelse(id %in% clusters[[1]],"Thalianol",
                                       ifelse(id %in% clusters[[2]],"Marneral",
                                              ifelse(id %in% clusters[[3]],"Tirucalladienol",
                                                     ifelse(id %in% clusters[[4]],"Arabidiol","None")))))

rownames(x) <- x$id
cluster_df = data.frame("Cluster" = x$cluster)
rownames(cluster_df) = rownames(x)
x <- select(x,-id,-cluster)
x <- log2(x+1)
Greens <- colorRampPalette(brewer.pal(n = 7, name ="Greens"))(100)
dendogram <- pheatmap(x,na_col = "black",cluster_rows = F,annotation_row = cluster_df,border_color = "NA",color = Greens)
ggsave(plot = dendogram,"~/iNOME/Figures/Main/Figure4/Parts/Dendogram_1_col.pdf")

