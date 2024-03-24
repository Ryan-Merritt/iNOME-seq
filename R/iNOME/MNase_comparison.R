rm(list=ls())
source("~/iNOME/Scripts/R/Libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
# Loading GFF -------------------------------------------------------------
gffs <- import("/Users/u1896702/work/ecotypes/annotations/blast_edit/TAIR10_genes_only.gff") %>%
  keepSeqlevels(.,paste0("Chr",seq(1,5,1)),"coarse")
# Loading RNA data --------------------------------------------------------
RNA <- read_delim("/Volumes/shared/Shared210/Ryan/iNOME/Ecotypes/ecotypes_rna/RNA/counts/edited/Col0_1_DEseq.input",delim = "\t",col_names = c("id","count"))
# Calculating TPM ---------------------------------------------------------
x <- gffs %>%
  filter(type == "exon") %>%
  mutate(id = unlist(gsub("\\..*","",Parent)))
lengths <- x %>%
  group_by(id) %>%
  group_split()
lengths <- lapply(lengths,function(x){
  x <- tibble(id = x$id[1],length = sum(width(reduce_ranges(x))))
  x
})
lengths <- bind_rows(lengths)

lengths <- lengths %>%
  inner_join(RNA) %>%
  mutate(TPM = count/length,TPM = TPM/(sum(TPM)/1000000))
# Loading TSS and ordering by expression ----------------------------------
tss <- gffs %>%
  mutate(start = ifelse(strand == "-",end,start),end = ifelse(strand == "+",start,end)) %>%
  filter(type == "gene")

tss <- as_tibble(tss) %>%
  distinct(ID,.keep_all = T) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

tss <- as_tibble(tss) %>%
  inner_join(lengths,by = c("ID" = "id")) %>%
  select(-length,-count) %>%
  arrange(desc(TPM)) %>%
  filter(TPM < mean(TPM))

tss <- tss %>%
  mutate(rank = seq(1,nrow(tss),1)) %>%
  makeGRangesFromDataFrame(keep.extra.columns = T)

tss_expressed <- tss %>%
  filter(TPM != 0) %>%
  mutate(quartile = ntile(TPM,3))
tss_zero <- tss %>%
  filter(TPM == 0) %>%
  mutate(quartile = 4)
tss <- bind_ranges(tss_expressed,tss_zero)
tss <- tss %>%
  mutate(quartile = ifelse(quartile == 1,"Low",
                           ifelse(quartile == 2,"Medium",
                                  ifelse(quartile == 3,"High",
                                         ifelse(quartile == 4,"Not Expressed",
                                                quartile)))))
names(tss) <- tss$rank
tss$quartile <- factor(tss$quartile,levels = c("High","Medium","Low","Not Expressed"))
tss <- dropSeqlevels(tss,c("ChrM","ChrC"),"coarse")

# Analysing data ----------------------------------------------------------
quartiles <- c("High","Medium","Low","Not Expressed")
mNASE_bw <- list()
iNOME_bw <- list()
ATAC_bw <- list()
for(i in quartiles){
  mNASE_bw[[i]] <- regionPlot(bamFile = "~/work/Julia/ecotypes_comparison/mnase/mNASE_col_50.bw",testRanges = filter(tss,quartile == i),format = "bigwig",paired = T,distanceUp = 500,distanceDown = 1000)
  iNOME_bw[[i]] <- regionPlot(bamFile = "~/work/Julia/ecotypes_comparison/moving_average/Col0_CHH_100.bw",testRanges = filter(tss,quartile == i),format = "bigwig",paired = T,distanceUp = 500,distanceDown = 1000)
  ATAC_bw[[i]] <- regionPlot(bamFile = "/Volumes/PhD_one/NOME_Seq/Col_features/bw/ATAC_col_10.bw",testRanges = filter(tss,quartile == i),format = "bigwig",paired = T,distanceUp = 500,distanceDown = 1000)
}
MNase_df <- list()
iNOME_df <- list()
ATAC_df <- list()
for(i in quartiles){
  MNase_df[[i]] <- assays(mNASE_bw[[i]])[[1]]
  iNOME_df[[i]] <- assays(iNOME_bw[[i]])[[1]]
  ATAC_df[[i]] <- assays(ATAC_bw[[i]])[[1]]
}

for(i in quartiles){
  MNase_df[[i]] <- tibble(Position = seq(1,1501,1),Score = colMeans(MNase_df[[i]]),Type = "MNase",quartile = i)
  iNOME_df[[i]] <- tibble(Position = seq(1,1501,1),Score = (1 - colMeans(iNOME_df[[i]])),Type = "iNOME",quartile = i)
  ATAC_df[[i]] <- tibble(Position = seq(1,1501,1),Score = (1 - colMeans(ATAC_df[[i]])),Type = "ATAC",quartile = i)
}

for(i in quartiles){
  MNase_df[[i]] <- tibble(Position = seq(1,1501,1),Score = colMeans(MNase_df[[i]]),Type = "MNase",quartile = i)
  iNOME_df[[i]] <- tibble(Position = seq(1,1501,1),Score = (1 - colMeans(iNOME_df[[i]]))-0.724,Type = "iNOME",quartile = i)
  ATAC_df[[i]] <- tibble(Position = seq(1,1501,1),Score = (1 - colMeans(ATAC_df[[i]]))-0.724,Type = "ATAC",quartile = i)
}
MNase_df <- bind_rows(MNase_df)
iNOME_df <- bind_rows(iNOME_df)
ATAC_df <- bind_rows(ATAC_df)
MNase_df <- MNase_df %>%
  mutate(norm = (Score-min(Score))/(max(Score)-min(Score)))
iNOME_df <- iNOME_df %>%
  mutate(norm = (Score-min(Score))/(max(Score)-min(Score)))
ATAC_df <- ATAC_df %>%
  mutate(norm = (Score-min(Score))/(max(Score)-min(Score)))

iNOME_MNase_comparison <- bind_rows(MNase_df,iNOME_df,ATAC_df)
iNOME_MNase_comparison$quartile <- factor(iNOME_MNase_comparison$quartile,levels = quartiles)

# Computing correlations --------------------------------------------------
results <- tibble(quartile = character(),p.value = numeric(),R = numeric())
for(i in quartiles){
  MNase <- iNOME_MNase_comparison %>%
    filter(quartile == i, Type == "MNase") %>%
    pull(norm)
  iNOME <- iNOME_MNase_comparison %>%
    filter(quartile == i, Type == "iNOME") %>%
    pull(norm)
  ATAC <- iNOME_MNase_comparison %>%
    filter(quartile == i, Type == "ATAC") %>%
    pull(norm)
  x <- cor.test(MNase,iNOME,method = "pearson")
  y <- cor.test(MNase,ATAC,method = "pearson")
  results <- bind_rows(results,tibble(quartile = i,p.value = x$p.value,R = x$estimate),tibble(quartile = i,p.value = y$p.value,R = y$estimate))
}
results <- results %>%
  mutate(x = 1200,y = -0.14,label = ifelse(p.value < 2.2e-16,paste0("p < ",2.2e-16),paste0("p = ",signif(p.value,2))),label = paste0("R = ",signif(R,2),", ",label),Method = rep(c("iNOME","ATAC"),4),significance = p.value<0.05)
# Plotting comparisons ----------------------------------------------------
iNOME_MNase_comparison <- iNOME_MNase_comparison %>%
  rename(Method = Type) %>%
  mutate(Method = ifelse(Method == "iNOME","1 - iNOME",ifelse(Method == "ATAC","1 - ATAC",Method)))
iNOME_MNase_comparison$quartile <- factor(iNOME_MNase_comparison$quartile,levels = c("High","Medium","Low","Not Expressed"))
x <- iNOME_MNase_comparison %>%
  #filter(Method != "1 - ATAC") %>%
ggplot() + 
  geom_line(aes(x = Position,y = norm,colour = Method)) + 
  geom_text(data = results %>% filter(Method == "iNOME"),aes(x = 1400,y = 0.025,label = label),show.legend = F) +
  geom_text(data = results %>% filter(Method == "ATAC"),aes(x = 1400,y = 0.15,label = label),show.legend = F) +
  facet_grid(rows = vars(quartile)) +
  scale_x_continuous(breaks = c(0,501,1001,1501),labels = c("-0.5kb","TSS","+0.5kb","+1kb")) +
  #scale_y_continuous("MNase",breaks = seq(0,0.2,0.05),labels = seq(0,0.2,0.05),sec.axis = sec_axis(~.+0.71,name = "1 - iNOME-seq",breaks = seq(0,0.9,0.1),labels = seq(0,0.9,0.1))) +
  scale_color_manual(values = c("#992025","#004fa7","#006500")) + 
  geom_vline(xintercept = 501,linetype = "dashed") +
  annotate("segment",x=-Inf, xend=Inf, y=-Inf, yend=-Inf,size = 1.2) +
  theme_cowplot()# +
  #theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title = element_blank())
ggsave(plot = x,"~/iNOME/Figures/Main/Figure3/Parts/test_2.pdf",height = 7,width  = 9)

