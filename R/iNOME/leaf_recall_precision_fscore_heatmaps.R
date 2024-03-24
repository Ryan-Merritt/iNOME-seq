rm(list=ls())
source("~/iNOME/Scripts/R/libraries.R")
source("~/iNOME/Scripts/R/Functions.R")
# Loading Leaf DMRs
files <- dir("/Volumes/PhD_one/NOME_Seq/data/DMRs/leaf/",pattern = "ME",full.names = T)
files <- files[grep("CM",files)]
for(i in files){
  load(i)
}

leaf_6h <- list(CG = leaf_6h_CpG_ME_CM,CHG = leaf_6h_CHG_ME_CM,CHH = leaf_6h_CHH_ME_CM)
for(i in 1:length(leaf_6h)){
  for(j in 1:length(leaf_6h[[i]])){
    for(k in 1:length(leaf_6h[[i]][[j]])){
      leaf_6h[[i]][[j]][[k]] <- leaf_6h[[i]][[j]][[k]] %>%
        {if(length(leaf_6h[[i]][[j]][[k]])>0) filter(.,regionType == "gain") else leaf_6h[[i]][[j]][[k]]}
    }
  }
}
leaf_12h <- list(CG = leaf_12h_CpG_ME_CM,CHG = leaf_12h_CHG_ME_CM,CHH = leaf_12h_CHH_ME_CM)
for(i in 1:length(leaf_12h)){
  for(j in 1:length(leaf_12h[[i]])){
    for(k in 1:length(leaf_12h[[i]][[j]])){
      leaf_12h[[i]][[j]][[k]] <- leaf_12h[[i]][[j]][[k]] %>%
      {if(length(leaf_12h[[i]][[j]][[k]])>0) filter(.,regionType == "gain") else leaf_12h[[i]][[j]][[k]]}
    }
  }
}
leaf_24h <- list(CG = leaf_24h_CpG_ME_CM,CHG = leaf_24h_CHG_ME_CM,CHH = leaf_24h_CHH_ME_CM)
for(i in 1:length(leaf_24h)){
  for(j in 1:length(leaf_24h[[i]])){
    for(k in 1:length(leaf_24h[[i]][[j]])){
      leaf_24h[[i]][[j]][[k]] <- leaf_24h[[i]][[j]][[k]] %>%
      {if(length(leaf_24h[[i]][[j]][[k]])>0) filter(.,regionType == "gain") else leaf_24h[[i]][[j]][[k]]}
    }
  }
}

rm(list = ls(pattern = "ME_CM"))
# Loading ATAC_seq data
ATAC_peaks <- import("/Volumes/PhD_one/NOME_Seq/data/features/ATAC_seq/seedling/GSM2260232_Col-0_50k_seedling_nuclei_FANS-ATAC-seq_rep2.regulatory_region.bed.gz")
seqlevels(ATAC_peaks) <- paste0("Chr",seq(1,5,1))

# Calculating TP, FP and FN
recall_6h <- vector("list",length(leaf_12h))
names(recall_6h) <- names(leaf_12h)
precision_6h <- vector("list",length(leaf_12h))
names(precision_6h) <- names(leaf_12h)
f_score_6h <- vector("list",length(leaf_12h))
names(f_score_6h) <- names(leaf_12h)
for(i in 1:length(leaf_6h)){
  recall_6h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_6h[[1]]))*100,"%"),names(leaf_6h[[1]][[1]])))
  precision_6h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_6h[[1]]))*100,"%"),names(leaf_6h[[1]][[1]])))
  f_score_6h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_6h[[1]]))*100,"%"),names(leaf_6h[[1]][[1]])))
  for(j in 1:length(leaf_6h[[i]])){
    for(k in 1:length(leaf_6h[[i]][[j]])){
      TP <- sum(width(intersect_ranges(leaf_6h[[i]][[j]][[k]],ATAC_peaks)))/sum(width(ATAC_peaks))
      FP <- sum(width(setdiff_ranges(leaf_6h[[i]][[j]][[k]],ATAC_peaks)))/sum(width(ATAC_peaks))
      FN <- sum(width(setdiff_ranges(ATAC_peaks,leaf_6h[[i]][[j]][[k]])))/sum(width(ATAC_peaks))
      recall_6h[[i]][j,k] <- (TP/(TP+FN))
      precision_6h[[i]][j,k] <- (TP/(TP+FP))
      f_score_6h[[i]][j,k] <- (2*TP/(2*TP+FP+FN))
    }
  }
  recall_6h[[i]][is.nan(recall_6h[[i]])] <- 0
  precision_6h[[i]][is.nan(precision_6h[[i]])] <- 0
  f_score_6h[[i]][is.nan(f_score_6h[[i]])] <- 0
  recall_6h[[i]] <- tibble(score = c(recall_6h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
         window = sort(as.numeric(rep(colnames(recall_6h[[i]]),10))),context = names(leaf_6h)[i],time = "6h")
  precision_6h[[i]] <- tibble(score = c(precision_6h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
                           window = sort(as.numeric(rep(colnames(precision_6h[[i]]),10))),context = names(leaf_6h)[i],time = "6h")
  f_score_6h[[i]] <- tibble(score = c(f_score_6h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
                           window = sort(as.numeric(rep(colnames(f_score_6h[[i]]),10))),context = names(leaf_6h)[i],time = "6h")
}
recall_6h <- bind_rows(recall_6h)
precision_6h <- bind_rows(precision_6h)
f_score_6h <- bind_rows(f_score_6h)

recall_12h <- vector("list",length(leaf_12h))
names(recall_12h) <- names(leaf_12h)
precision_12h <- vector("list",length(leaf_12h))
names(precision_12h) <- names(leaf_12h)
f_score_12h <- vector("list",length(leaf_12h))
names(f_score_12h) <- names(leaf_12h)
for(i in 1:length(leaf_12h)){
  recall_12h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_12h[[1]]))*100,"%"),names(leaf_12h[[1]][[1]])))
  precision_12h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_12h[[1]]))*100,"%"),names(leaf_12h[[1]][[1]])))
  f_score_12h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_12h[[1]]))*100,"%"),names(leaf_12h[[1]][[1]])))
  for(j in 1:length(leaf_12h[[i]])){
    for(k in 1:length(leaf_12h[[i]][[j]])){
      TP <- sum(width(intersect_ranges(leaf_12h[[i]][[j]][[k]],ATAC_peaks)))/sum(width(ATAC_peaks))
      FP <- sum(width(setdiff_ranges(leaf_12h[[i]][[j]][[k]],ATAC_peaks)))/sum(width(ATAC_peaks))
      FN <- sum(width(setdiff_ranges(ATAC_peaks,leaf_12h[[i]][[j]][[k]])))/sum(width(ATAC_peaks))
      recall_12h[[i]][j,k] <- (TP/(TP+FN))
      precision_12h[[i]][j,k] <- (TP/(TP+FP))
      f_score_12h[[i]][j,k] <- (2*TP/(2*TP+FP+FN))
    }
  }
  recall_12h[[i]][is.nan(recall_12h[[i]])] <- 0
  precision_12h[[i]][is.nan(precision_12h[[i]])] <- 0
  f_score_12h[[i]][is.nan(f_score_12h[[i]])] <- 0
  recall_12h[[i]] <- tibble(score = c(recall_12h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
                            window = sort(as.numeric(rep(colnames(recall_12h[[i]]),10))),context = names(leaf_12h)[i],time = "12h")
  precision_12h[[i]] <- tibble(score = c(precision_12h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
                               window = sort(as.numeric(rep(colnames(precision_12h[[i]]),10))),context = names(leaf_12h)[i],time = "12h")
  f_score_12h[[i]] <- tibble(score = c(f_score_12h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
                             window = sort(as.numeric(rep(colnames(f_score_12h[[i]]),10))),context = names(leaf_12h)[i],time = "12h")
}
recall_12h <- bind_rows(recall_12h)
precision_12h <- bind_rows(precision_12h)
f_score_12h <- bind_rows(f_score_12h)

recall_24h <- vector("list",length(leaf_24h))
names(recall_24h) <- names(leaf_24h)
precision_24h <- vector("list",length(leaf_24h))
names(precision_24h) <- names(leaf_24h)
f_score_24h <- vector("list",length(leaf_24h))
names(f_score_24h) <- names(leaf_24h)
for(i in 1:length(leaf_24h)){
  recall_24h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_24h[[1]]))*100,"%"),names(leaf_24h[[1]][[1]])))
  precision_24h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_24h[[1]]))*100,"%"),names(leaf_24h[[1]][[1]])))
  f_score_24h[[i]] <- matrix(ncol = 15,nrow = 10,dimnames = list(paste0(as.numeric(names(leaf_24h[[1]]))*100,"%"),names(leaf_24h[[1]][[1]])))
  for(j in 1:length(leaf_24h[[i]])){
    for(k in 1:length(leaf_24h[[i]][[j]])){
      TP <- sum(width(intersect_ranges(leaf_24h[[i]][[j]][[k]],ATAC_peaks)))/sum(width(ATAC_peaks))
      FP <- sum(width(setdiff_ranges(leaf_24h[[i]][[j]][[k]],ATAC_peaks)))/sum(width(ATAC_peaks))
      FN <- sum(width(setdiff_ranges(ATAC_peaks,leaf_24h[[i]][[j]][[k]])))/sum(width(ATAC_peaks))
      recall_24h[[i]][j,k] <- (TP/(TP+FN))
      precision_24h[[i]][j,k] <- (TP/(TP+FP))
      f_score_24h[[i]][j,k] <- (2*TP/(2*TP+FP+FN))
    }
  }
  recall_24h[[i]][is.nan(recall_24h[[i]])] <- 0
  precision_24h[[i]][is.nan(precision_24h[[i]])] <- 0
  f_score_24h[[i]][is.nan(f_score_24h[[i]])] <- 0
  recall_24h[[i]] <- tibble(score = c(recall_24h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
                           window = sort(as.numeric(rep(colnames(recall_24h[[i]]),10))),context = names(leaf_24h)[i],time = "24h")
  precision_24h[[i]] <- tibble(score = c(precision_24h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
                              window = sort(as.numeric(rep(colnames(precision_24h[[i]]),10))),context = names(leaf_24h)[i],time = "24h")
  f_score_24h[[i]] <- tibble(score = c(f_score_24h[[i]]),methylation_difference = as.character(rep(c(10,15,seq(20,90,10)),15)),
                            window = sort(as.numeric(rep(colnames(f_score_24h[[i]]),10))),context = names(leaf_24h)[i],time = "24h")
}
recall_24h <- bind_rows(recall_24h)
precision_24h <- bind_rows(precision_24h)
f_score_24h <- bind_rows(f_score_24h)

recall <- bind_rows(recall_6h,recall_12h,recall_24h)
precision <- bind_rows(precision_6h,precision_12h,precision_24h)
f_score <- bind_rows(f_score_6h,f_score_12h,f_score_24h)

recall$methylation_difference <- factor(recall$methylation_difference,levels = c("10", "15", "20", "30", "40", "50", "60", "70", "80", "90"))
recall$window <- factor(recall$window,levels = c("50", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1500", "2000", "3000", "5000"))
recall$time <- factor(recall$time,levels = c("6h","12h","24h"))
recall <- recall %>%
  mutate(label = round(score,1))

precision$methylation_difference <- factor(precision$methylation_difference,levels = c("10", "15", "20", "30", "40", "50", "60", "70", "80", "90"))
precision$window <- factor(precision$window,levels = c("50", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1500", "2000", "3000", "5000"))
precision$time <- factor(precision$time,levels = c("6h","12h","24h"))
precision <- precision %>%
  mutate(label = round(score,1))

f_score$methylation_difference <- factor(f_score$methylation_difference,levels = c("10", "15", "20", "30", "40", "50", "60", "70", "80", "90"))
f_score$window <- factor(f_score$window,levels = c("50", "100", "200", "300", "400", "500", "600", "700", "800", "900", "1000", "1500", "2000", "3000", "5000"))
f_score$time <- factor(f_score$time,levels = c("6h","12h","24h"))
f_score <- f_score %>%
  mutate(label = round(score,1))
recall %>%
  mutate(label = ifelse(label == 0,NA,label)) %>%
ggplot(aes(x = window,y = methylation_difference,fill = score)) +
  geom_tile() +
  scale_fill_gradient(low = '#ffffff',high = '#006cff') +
  facet_grid(cols = vars(context),rows = vars(time)) +
  geom_text(aes(label = label)) +
  annotate("segment",x = -Inf,xend = Inf,y = -Inf,yend = -Inf,size = 1) +
  annotate("segment",x = -Inf,xend = Inf,y = Inf,yend = Inf,size = 1) +
  annotate("segment",x = -Inf,xend = -Inf,y = -Inf,yend = Inf,size = 1) +
  annotate("segment",x = Inf,xend = Inf,y = -Inf,yend = Inf,size = 1) +
  xlab("Window") +
  ylab("Methylation Difference") +
  theme_cowplot() +
  theme(strip.background = element_blank(),strip.text.y = element_text(angle = 0),axis.text.x = element_text(angle = 90,vjust = 0.5))
ggsave("~/iNOME/Figures/Supplementary/leaf_recall_heatmap_no_0.pdf",width = 15)

precision %>%
  mutate(label = ifelse(label == 0,NA,label)) %>%
ggplot(aes(x = window,y = methylation_difference,fill = score)) +
  geom_tile() +
  scale_fill_gradient(low = '#ffffff',high = '#006cff') +
  facet_grid(cols = vars(context),rows = vars(time)) +
  geom_text(aes(label = label)) +
  annotate("segment",x = -Inf,xend = Inf,y = -Inf,yend = -Inf,size = 1) +
  annotate("segment",x = -Inf,xend = Inf,y = Inf,yend = Inf,size = 1) +
  annotate("segment",x = -Inf,xend = -Inf,y = -Inf,yend = Inf,size = 1) +
  annotate("segment",x = Inf,xend = Inf,y = -Inf,yend = Inf,size = 1) +
  xlab("Window") +
  ylab("Methylation Difference") +
  theme_cowplot() +
  theme(strip.background = element_blank(),strip.text.y = element_text(angle = 0),axis.text.x = element_text(angle = 90,vjust = 0.5))
ggsave("~/iNOME/Figures/Supplementary/leaf_precision_heatmap_no_0.pdf",width = 15)

f_score %>%
  mutate(label = ifelse(label == 0,NA,label)) %>%
ggplot(aes(x = window,y = methylation_difference,fill = score)) +
  geom_tile() +
  scale_fill_gradient(low = '#ffffff',high = '#006cff') +
  facet_grid(cols = vars(context),rows = vars(time)) +
  geom_text(aes(label = label)) +
  annotate("segment",x = -Inf,xend = Inf,y = -Inf,yend = -Inf,size = 1) +
  annotate("segment",x = -Inf,xend = Inf,y = Inf,yend = Inf,size = 1) +
  annotate("segment",x = -Inf,xend = -Inf,y = -Inf,yend = Inf,size = 1) +
  annotate("segment",x = Inf,xend = Inf,y = -Inf,yend = Inf,size = 1) +
  xlab("Window") +
  ylab("Methylation Difference") +
  theme_cowplot() +
  theme(strip.background = element_blank(),strip.text.y = element_text(angle = 0),axis.text.x = element_text(angle = 90,vjust = 0.5))
ggsave("~/iNOME/Figures/Supplementary/leaf_f_score_heatmap_no_0.pdf",width = 15)

