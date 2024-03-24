source("~/iNOME/Scripts/R/libraries.R")
# Root leaf annotation ----------------------------------------------------
gff <- import("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff") %>%
  dropSeqlevels(.,c("ChrC","ChrM"),"coarse")
# Extracting genes
genes <- filter(gff,type == "gene")
# Extracting TEs that dont overlap with genes
TE <- filter(gff,type == "transposable_element")
TE <- setdiff_ranges(TE,genes)
TE <- reduce_ranges(TE)
# Extracting regions upstream of genes that dont overlap with other features
one_kb_upstream <- setdiff_ranges(promoters(genes,upstream = 1000,downstream = 0),genes)
one_kb_upstream <- setdiff_ranges(one_kb_upstream,TE)
# Extracting regions downstream of genes that dont overlap with other features
one_kb_downstream <- genes
start(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"])
end(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"]) + 1000
end(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"])
start(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"]) - 1000

one_kb_downstream <- setdiff_ranges(one_kb_downstream,genes)
one_kb_downstream <- setdiff_ranges(one_kb_downstream,TE)
one_kb_downstream <- setdiff_ranges(one_kb_downstream,one_kb_upstream)
# checking removal was correct
one_kb_upstream <- setdiff_ranges(one_kb_upstream,one_kb_downstream)
one_kb_downstream <- setdiff_ranges(one_kb_downstream,one_kb_upstream)
# Calculating widths of overlaps
## Total expected
widths <- c(sum(width(reduce_ranges(genes))),sum(width(reduce_ranges(TE))),sum(width(reduce_ranges(one_kb_upstream))),sum(width(reduce_ranges(one_kb_downstream))))
total_overlap <- tibble(Type = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"),
                            Length = c(widths,sum(width(getWholeChromosomes(gff))) - sum(widths)),
                            Sample = "Expected") %>%
  mutate(Percentage = (Length/sum(Length))*100)
## Root overlap
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/root/CHH_ME_CE_12h.RData")

DMRs_root <- root_12h_CHH_ME_CE$`0.3`$`700`

root_widths <- c()
root_widths <- c(root_widths,sum(width(reduce_ranges(intersect_ranges(DMRs_root,genes)))))
DMRs_root <- setdiff_ranges(DMRs_root,genes)
root_widths <- c(root_widths,sum(width(reduce_ranges(intersect_ranges(DMRs_root,TE)))))
DMRs_root <- setdiff_ranges(DMRs_root,TE)
root_widths <- c(root_widths,sum(width(reduce_ranges(intersect_ranges(DMRs_root,one_kb_upstream)))))
DMRs_root <- setdiff_ranges(DMRs_root,one_kb_upstream)
root_widths <- c(root_widths,sum(width(reduce_ranges(intersect_ranges(DMRs_root,one_kb_downstream)))))
DMRs_root <- setdiff_ranges(DMRs_root,one_kb_downstream)
root_widths <- c(root_widths,sum(width(reduce_ranges(DMRs_root))))
root_overlap <- tibble(Type = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"),
                           Length = c(root_widths),
                           Sample = "Root") %>%
  mutate(Percentage = (Length/sum(Length))*100)
#sum(width(getWholeChromosomes(gff))) - sum(root_widths)
## Leaf overlap
load("/Volumes/PhD_one/NOME_Seq/data/DMRs/leaf/CHH_ME_CE_12h.RData")

DMRs_leaf <- leaf_12h_CHH_ME_CE$`0.15`$`700`

leaf_widths <- c()
leaf_widths <- c(leaf_widths,sum(width(reduce_ranges(intersect_ranges(DMRs_leaf,genes)))))
DMRs_leaf <- setdiff_ranges(DMRs_leaf,genes)
leaf_widths <- c(leaf_widths,sum(width(reduce_ranges(intersect_ranges(DMRs_leaf,TE)))))
DMRs_leaf <- setdiff_ranges(DMRs_leaf,TE)
leaf_widths <- c(leaf_widths,sum(width(reduce_ranges(intersect_ranges(DMRs_leaf,one_kb_upstream)))))
DMRs_leaf <- setdiff_ranges(DMRs_leaf,one_kb_upstream)
leaf_widths <- c(leaf_widths,sum(width(reduce_ranges(intersect_ranges(DMRs_leaf,one_kb_downstream)))))
DMRs_leaf <- setdiff_ranges(DMRs_leaf,one_kb_downstream)
leaf_widths <- c(leaf_widths,sum(width(reduce_ranges(DMRs_leaf))))
leaf_overlap <- tibble(Type = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"),
                           Length = c(leaf_widths),
                           Sample = "Leaf") %>%
  mutate(Percentage = (Length/sum(Length))*100)
# Plotting graph
all_overlaps <- bind_rows(total_overlap,root_overlap,leaf_overlap) %>%
  group_by(Sample) %>%
  mutate(total = sum(Length),Percentage = (Length/total)*100)

all_overlaps$Type <- factor(all_overlaps$Type,levels = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"))

all_overlaps$Sample <- factor(all_overlaps$Sample,levels = c("Root", "Leaf", "Expected"))

# '#005800', '#038202', '#45aa34', '#71d45c', '#a0ff87' green palatte
# '#3a218a', '#645186', '#808080', '#965140', '#980000' Blue red
# '#3a218a', '#645186', '#808080', '#608b54', '#008000' blue green
# c('#980000', '#965140', '#808080', '#a37bbe', '#401f59') Red purple
(p1 <- ggplot(all_overlaps,aes(x = Sample,y = Percentage,fill = Type)) +
  geom_bar(stat = "identity",colour = "black",size = 0.1) +
  scale_fill_manual(values = c('#3a218a', '#645186', '#808080', '#965140', '#980000')) +
  theme_cowplot() +
    theme(axis.line = element_line(size = 0.1)))
ggsave("~/iNOME/Figures/Main/Figure2/Parts/root_length_distribution_nooverlap.pdf",width = 10,height = 6)
# calculating log2(observed/expected)
#colours <- colorRampPalette(c("#3a218a","white","#980000"))(5)
obs_exp_tab <- tibble(Type = rep(root_overlap$Type,2),Observed = c(root_overlap$Percentage,leaf_overlap$Percentage),
                      Expected = rep(total_overlap$Percentage,2),Sample = c(rep("Root",5),rep("Leaf",5))) %>%
  mutate(Stat = log2(Observed/Expected))

split_tab <- obs_exp_tab %>%
  group_split(Sample)
split_tab <- tibble(Type = split_tab[[1]]$Type,Root = split_tab[[2]]$Stat,Leaf = split_tab[[1]]$Stat) %>%
  as.data.frame()
rownames(split_tab) <- split_tab$Type
split_tab <- select(split_tab,-Type)
split_tab <- t(split_tab)
heatmap <- pheatmap(split_tab,cluster_rows = F,cluster_cols = F,border_color = "NA",display_numbers = T,
                    cellwidth = 50,cellheight = 50)
ggsave(plot = heatmap,"./iNOME/Figures/Main/Figure2/Parts/root_length_log2.pdf")
# common and unique DMRs --------------------------------------------------
gff <- import("~/work/ecotypes/annotations/Col0/TAIR10_genes_transposons.gff") %>%
  dropSeqlevels(.,c("ChrC","ChrM"),"coarse")
# Extracting genes
genes <- filter(gff,type == "gene")
# Extracting TEs that dont overlap with genes
TE <- filter(gff,type == "transposable_element")
TE <- setdiff_ranges(TE,genes)
TE <- reduce_ranges(TE)
# Extracting regions upstream of genes that dont overlap with other features
one_kb_upstream <- setdiff_ranges(promoters(genes,upstream = 1000,downstream = 0),genes)
one_kb_upstream <- setdiff_ranges(one_kb_upstream,TE)
# Extracting regions downstream of genes that dont overlap with other features
one_kb_downstream <- genes
start(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"])
end(one_kb_downstream[strand(one_kb_downstream)=="+"]) <- end(one_kb_downstream[strand(one_kb_downstream)=="+"]) + 1000
end(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"])
start(one_kb_downstream[strand(one_kb_downstream)=="-"]) <- start(one_kb_downstream[strand(one_kb_downstream)=="-"]) - 1000

one_kb_downstream <- setdiff_ranges(one_kb_downstream,genes)
one_kb_downstream <- setdiff_ranges(one_kb_downstream,TE)
one_kb_downstream <- setdiff_ranges(one_kb_downstream,one_kb_upstream)
# checking removal was correct
one_kb_upstream <- setdiff_ranges(one_kb_upstream,one_kb_downstream)
one_kb_downstream <- setdiff_ranges(one_kb_downstream,one_kb_upstream)
# Calculating widths of overlaps
## Total expected
widths <- c(sum(width(reduce_ranges(genes))),sum(width(reduce_ranges(TE))),sum(width(reduce_ranges(one_kb_upstream))),sum(width(reduce_ranges(one_kb_downstream))))
total_overlap <- tibble(Type = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"),
                            Length = c(widths,sum(width(getWholeChromosomes(gff))) - sum(widths)),
                            Sample = "Expected") %>%
  mutate(Percentage = (Length/sum(Length))*100)

# Loading ATAC-seq peaks
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/Col_features/peaks/ATAC_col_q005_peaks.narrowPeak")
ATAC_peaks <- keepSeqlevels(ATAC_peaks,c("Chr1", "Chr2", "Chr3", "Chr4", "Chr5"),"coarse")
ATAC_peaks <- GRanges(seqnames(ATAC_peaks),IRanges(start(ATAC_peaks),end(ATAC_peaks)))
ATAC_peaks <- unique(ATAC_peaks)
# Load DMRs
DMRs <- root_12h_CHH_ME_CE$`0.3`$`700`
cDMRs <- intersect_ranges(DMRs,ATAC_peaks)

## Extracting unique DMRs
sDMRs <- setdiff_ranges(DMRs,ATAC_peaks)
sATAC <- setdiff_ranges(ATAC_peaks,DMRs)

## cDMRs overlap
cDMRs_widths <- c()
cDMRs_widths <- c(cDMRs_widths,sum(width(reduce_ranges(intersect_ranges(cDMRs,genes)))))
cDMRs <- setdiff_ranges(cDMRs,genes)
cDMRs_widths <- c(cDMRs_widths,sum(width(reduce_ranges(intersect_ranges(cDMRs,TE)))))
cDMRs <- setdiff_ranges(cDMRs,TE)
cDMRs_widths <- c(cDMRs_widths,sum(width(reduce_ranges(intersect_ranges(cDMRs,one_kb_upstream)))))
cDMRs <- setdiff_ranges(cDMRs,one_kb_upstream)
cDMRs_widths <- c(cDMRs_widths,sum(width(reduce_ranges(intersect_ranges(cDMRs,one_kb_downstream)))))
cDMRs <- setdiff_ranges(cDMRs,one_kb_downstream)
cDMRs_widths <- c(cDMRs_widths,sum(width(reduce_ranges(cDMRs))))
cDMRs_overlap <- tibble(Type = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"),
                            Length = c(cDMRs_widths),
                            Sample = "Common") %>%
  mutate(Percentage = (Length/sum(Length)*100))
#sum(width(getWholeChromosomes(gff))) - sum(cDMRs_widths)
## sDMRs overlap
sDMRs_widths <- c()
sDMRs_widths <- c(sDMRs_widths,sum(width(reduce_ranges(intersect_ranges(sDMRs,genes)))))
sDMRs <- setdiff_ranges(sDMRs,genes)
sDMRs_widths <- c(sDMRs_widths,sum(width(reduce_ranges(intersect_ranges(sDMRs,TE)))))
sDMRs <- setdiff_ranges(sDMRs,TE)
sDMRs_widths <- c(sDMRs_widths,sum(width(reduce_ranges(intersect_ranges(sDMRs,one_kb_upstream)))))
sDMRs <- setdiff_ranges(sDMRs,one_kb_upstream)
sDMRs_widths <- c(sDMRs_widths,sum(width(reduce_ranges(intersect_ranges(sDMRs,one_kb_downstream)))))
sDMRs <- setdiff_ranges(sDMRs,one_kb_downstream)
sDMRs_widths <- c(sDMRs_widths,sum(width(reduce_ranges(sDMRs))))
sDMRs_overlap <- tibble(Type = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"),
                            Length = c(sDMRs_widths),
                            Sample = "iNOME") %>%
  mutate(Percentage = (Length/sum(Length)*100))
# ATAC overlap
sATAC_widths <- c()
sATAC_widths <- c(sATAC_widths,sum(width(reduce_ranges(intersect_ranges(sATAC,genes)))))
sATAC <- setdiff_ranges(sATAC,genes)
sATAC_widths <- c(sATAC_widths,sum(width(reduce_ranges(intersect_ranges(sATAC,TE)))))
sATAC <- setdiff_ranges(sATAC,TE)
sATAC_widths <- c(sATAC_widths,sum(width(reduce_ranges(intersect_ranges(sATAC,one_kb_upstream)))))
sATAC <- setdiff_ranges(sATAC,one_kb_upstream)
sATAC_widths <- c(sATAC_widths,sum(width(reduce_ranges(intersect_ranges(sATAC,one_kb_downstream)))))
sATAC <- setdiff_ranges(sATAC,one_kb_downstream)
sATAC_widths <- c(sATAC_widths,sum(width(reduce_ranges(sATAC))))
sATAC_overlap <- tibble(Type = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"),
                        Length = c(sATAC_widths),
                        Sample = "ATAC") %>%
  mutate(Percentage = (Length/sum(Length)*100))
# Plotting graph
all_overlaps <- bind_rows(total_overlap,cDMRs_overlap,sDMRs_overlap,sATAC_overlap) %>%
  group_by(Sample) %>%
  mutate(total = sum(Length),Percentage = (Length/total)*100)

all_overlaps$Type <- factor(all_overlaps$Type,levels = c("Genes","Transposable Elements","1kb Upstream","1kb Downstream","Intergenic"))

all_overlaps$Sample <- factor(all_overlaps$Sample,levels = c("Common", "iNOME","ATAC", "Expected"))

(p1 <- ggplot(all_overlaps,aes(x = Sample,y = Percentage,fill = Type)) +
  geom_bar(stat = "identity",colour = "black",size = 0.1) +
  scale_fill_manual(values = c('#3a218a', '#645186', '#808080', '#965140', '#980000')) +
  theme_cowplot() +
  theme(axis.line = element_line(size = 0.1)))
ggsave("~/iNOME/Figures/Main/Figure2/Parts/common_unique_distribution_nooverlap.pdf",width = 10,height = 6)
# calculating log2(observed/expected)
#'#3a218a', '#645186', '#808080', '#965140', '#980000'
colours <- colorRampPalette(c("#3B4992FF","white","#DC0000FF"))(10)
obs_exp_tab <- tibble(Type = rep(cDMRs_overlap$Type,3),Observed = c(cDMRs_overlap$Percentage,sDMRs_overlap$Percentage,sATAC_overlap$Percentage),
                      Expected = rep(total_overlap$Percentage,3),Sample = c(rep("Common",5),rep("iNOME",5),rep("ATAC",5))) %>%
  mutate(Stat = log2(Observed/Expected))
obs_exp_tab$Sample <- factor(obs_exp_tab$Sample,levels = c("ATAC","iNOME","Common"))
ggplot(obs_exp_tab,aes(x = Type,y = Sample,fill = Stat)) +
  geom_tile(colour = "black",size = 0.1) +
  coord_equal() +
  scale_fill_gradientn(colours = colours,limits = c(-2,1.5)) +
  theme_cowplot() +
  theme(axis.title = element_blank(),axis.ticks = element_blank(),axis.line = element_blank())
ggsave("~/iNOME/Figures/Main/Figure2/Parts/common_unique_log2.pdf")
