# this script analyses the TSL read count data obtained from featureCounts
# the inputs are the RDS DESEQ2 data objects generated by the script gene_annotation_correction.R

#
# load RDS objects (from many different datasets)
#
dsets <- c("PRJNA12603.featureCounts.10bp.gene.RDS",
           "STRINGTIE.featureCounts.10bp.gene.RDS",
           "BRAKER_BAM.featureCounts.10bp.gene.RDS",
           "BRAKER_TRINITY.featureCounts.10bp.gene.RDS",
           "DENOVO_ALL.featureCounts.10bp.gene.RDS",
           "PRJNA12603.featureCounts.8bp.gene.RDS",
           "STRINGTIE.featureCounts.8bp.gene.RDS",
           "BRAKER_BAM.featureCounts.8bp.gene.RDS",
           "BRAKER_TRINITY.featureCounts.8bp.gene.RDS",
           "DENOVO_ALL.featureCounts.8bp.gene.RDS")
dd.genes <- lapply(paste("datasets/", dsets, sep=""), readRDS)
names(dd.genes) <- dsets

dsets <- c("PRJNA12603.featureCounts.10bp.exon.conv.RDS",
           "STRINGTIE.featureCounts.10bp.exon.conv.RDS",
           "BRAKER_BAM.featureCounts.10bp.exon.conv.RDS",
           "BRAKER_TRINITY.featureCounts.10bp.exon.conv.RDS",
           "DENOVO_ALL.featureCounts.10bp.exon.conv.RDS",
           "PRJNA12603.featureCounts.8bp.exon.conv.RDS",
           "STRINGTIE.featureCounts.8bp.exon.conv.RDS",
           "BRAKER_BAM.featureCounts.8bp.exon.conv.RDS",
           "BRAKER_TRINITY.featureCounts.8bp.exon.conv.RDS",
           "DENOVO_ALL.featureCounts.8bp.exon.conv.RDS")
dd.exons <- lapply(paste("datasets/", dsets, sep=""), readRDS)
names(dd.exons) <- dsets

#
# Multivariate exploration using Jaccard and normalized counts
#

dendrogram <- function(dis, an, xl=12){
  hc <- hclust(dis, method="ward.D2")
  ph <- as.phylo(hc)
  tiplabs <- ph$tip.label
  try(ph <- phytools::reroot(ph, MRCA(ph, ph$tip.label[grep("TSL2|TSL10|TSL12", ph$tip.label)]), position=0.5), silent=T)
  ph <- groupOTU(ph, split(tiplabs, an$Group2), "tsl")
  
  attr(ph, "tsl")[which(attr(ph, "tsl")==0)] <- "others"
  attr(ph, "tsl") <- factor(attr(ph, "tsl"))
  ph$tip.label <- gsub("LIB.*_", "", ph$tip.label)
  an$Label <- gsub(".*TSL", "Tsp-SL", an$ID)
  an$Label[an$Label=="Tsp-SL1345"] <- "Tsp-SL13+"

  tr <- ggtree(ph, aes(color=tsl), layout="rectangular", branch.length="none") +
    geom_tiplab(aes(label=Label), offset=0.1, size=3) + xlim(c(0,10))

  tr <- tr %<+% as.data.frame(an[,-1]) + 
    geom_tippoint(aes(shape=Group, fill=tsl), size=3) +
    scale_fill_manual(values = unname(pcols.tsl[c("OSL", "ISL")])) +
    scale_colour_manual(values = unname(pcols.tsl[c("OSL", "ISL")])) +
    scale_shape_manual(values=c(others=21, TSL2=22, TSL10=25, TSL12=24))
  return(tr)
}

ordinate <- function(vc, an, ll, tt){
  pdata <- data.frame(X1=vc[,1], X2=vc[,2], an)
  pdata$Group <- factor(pdata$Group, levels=levels(pdata$Group)[c(4,2,3,1)])
  lbl <- c(expression(paste(italic(" Tsp-"), "SL2  ")), 
           expression(paste(italic(" Tsp-"), "SL10  ")), 
           expression(paste(italic(" Tsp-"), "SL12  ")), "others")
  ord <- ggplot(pdata, aes(x=X1, y=X2, fill=Group, shape=Group)) + 
    geom_hline(yintercept=0, linetype="dashed", colour="grey") +
    geom_vline(xintercept=0, linetype="dashed", colour="grey") +
    geom_point(colour="white", size=5) +
    scale_fill_manual(values=c(others="dodgerblue3", TSL2="#FF9900", TSL10="#FF9900", TSL12="#FF9900"),
                      labels=lbl) + 
    scale_shape_manual(values=c(others=21, TSL2=22, TSL10=25, TSL12=24), 
                       labels=lbl) + 
    xlab(paste(ll, "1", sep="")) + ylab(paste(ll, "2", sep="")) + labs(title=tt) +
    #coord_equal() +
    theme_bw() +
    theme(legend.position="bottom", 
          legend.text = element_text(size=11),
          #legend.title = element_text(size=11),
          legend.title = element_blank())
  return(ord)
}

Jaccard <- function(des){
  des$TSL <- factor(des$TSL)
  des$LIB <- factor(des$LIB)
  
  # Jaccard distance 
  prab <- t(counts(des, normalized=F))
  prab[prab>0] <- 1
  jac <- dist.binary(prab, method=1)
  
  # MDS plot
  cm <- cmdscale(jac)
  mds <- ordinate(cm, colData(des), "MDS", "Jaccard distance")
  
  # dendrogram
  jac.tr <- dendrogram(jac, colData(des), 9.7)
  
  # DAPC
  dapc.plots <- DAPC(prab, colData(des)$Group2)
  
  return(list(mds, jac.tr, dapc.plots[[1]], dapc.plots[[2]], dapc.plots[[3]]))
}

normCounts <- function(des){
  # normalize counts
  des <- estimateSizeFactors(des)
  vsd <- varianceStabilizingTransformation(des, blind = TRUE, fitType = "parametric")
  vsdd <- dist(t(assay(vsd)))

  # PCA
  pc <- prcomp(t(assay(vsd)))
  pca <- ordinate(pc$x, colData(des), "PC", "Normalized-counts")
  
  # dendrogram
  vsdd.tr <- dendrogram(vsdd, colData(des), xl=13)
  
  # DAPC
  dapc.plots <- DAPC(t(assay(vsd)), colData(des)$Group2)
  
  return(list(pca, vsdd.tr, dapc.plots[[1]], dapc.plots[[2]], dapc.plots[[3]]))
}

DAPC <- function(da, gr){
  pc.best <- xvalDapc(da, gr, parallel = "snow", ncpus=4, xval.plot=F, scale=F)$`Number of PCs Achieving Lowest MSE`
  #pc.best=15
  dp <- list()
  dp$DAPC <- dapc(da, gr, n.pca=as.numeric(pc.best), n.da=1, var.loadings=T, scale=F, center=T)
  
  # density plot
  pdata <- data.frame(dp$DAPC$ind.coord, Group=gr)
  pdata$Group <- factor(pdata$Group, levels=rev(levels(pdata$Group)))
  dpp <- ggplot(pdata, aes(x=LD1, colour=Group, fill=Group)) + 
    geom_density(colour="white") +
    geom_rug() +
    geom_vline(xintercept=0, linetype="dashed", colour="grey") +
    #scale_x_continuous(expand=c(0,0)) +
    scale_fill_manual(values=unname(pcols.tsl[c("ISL", "OSL")]), guide=FALSE) +
    scale_color_manual(values=unname(pcols.tsl[c("ISL", "OSL")]), guide=FALSE) +
    xlab("LD1") + ylab("Density") +
    theme_bw()
  
  #plot(dp$DAPC$var.load, dp$DAPC$var.contr)
  
  # loading plot
  km <- kmeans(as.numeric(dp$DAPC$var.contr), 2)
  pdata <- data.frame(x=1:length(dp$DAPC$var.load), 
                      dp$DAPC$var.load, Cluster=km$cluster, 
                      ifelse(dp$DAPC$var.load>0, "ISL", "OSL"))
  #pdata$LD1.1[pdata$Cluster==2] <- NA
  lp <- ggplot(pdata, aes(x=x, xend=x, y=0, yend=LD1, colour=LD1.1)) + 
    geom_segment() +
    geom_hline(yintercept=0, linetype="dashed", colour="grey") +
    #scale_colour_gradientn(colours=c("lightgray", "black")[order(km$centers)], guide=FALSE) +
    scale_colour_manual(values=pcols.tsl, guide=FALSE) +
    xlab("Gene") + ylab("Contribution") +
    theme_bw()
  
  return(list(dpp, lp, pdata))
}

# run all multivariate analyses and make compound figure
multivar_figure <- function(des, fig){
  plots1 <- Jaccard(des)
  plots2 <- normCounts(des)
  
  # main figure is combined with Jonathan's RNA loop figure
  img <- png::readPNG("../analysis2/Marius_et_al_Figure_1_2.png")
  g <- grid::rasterGrob(img, interpolate=TRUE)
  #grid.draw(g + textGrob("A",gp=gpar(fontsize=20,font=2), x=0, y=0, hjust=0))
  
  fig.main <- grid.arrange(
               arrangeGrob(top=textGrob("A",gp=gpar(fontsize=20,font=2), x=0, hjust=0), 
                           g),  
               arrangeGrob(top=textGrob("B",gp=gpar(fontsize=20,font=2), x=0, hjust=0), 
                           plots1[[1]], padding = unit(0, "line")), 
               arrangeGrob(top=textGrob("C",gp=gpar(fontsize=20,font=2), x=0, hjust=0), 
                           plots1[[2]], padding = unit(0, "line"), clip="on"), 
               arrangeGrob(top=textGrob("D",gp=gpar(fontsize=20,font=2), x=0, hjust=0),
                           plots1[[3]], plots1[[4]], nrow=2, clip="on"),
               layout_matrix=matrix(c(1,1,1,2,3,4), ncol=3, byrow=T))
  
  
  
  # using multipanelfigure
  #fig.main <- multi_panel_figure(width=c(4,4,4), height=c(3,4), unit="inches",
  #                              column_spacing = 0, row_spacing = 0.3)
                                
  #fig.main %<>% fill_panel("../analysis2/Marius_et_al_Figure_1_2.png", row = 1, column=c(1,2,3), scaling="fit")
  #fig.main %<>% fill_panel(plots1[[1]], row = 1, column=1)
  #fig.main %<>% fill_panel(plots1[[2]], row = 1, column=2)
  #fig.main %<>% fill_panel(arrangeGrob(plots1[[3]], plots1[[4]], nrow=2, clip="on"), row = 1, column=3)
  
  #fig.main %<>% fill_panel(plots1[[1]], row = 2, column=1, scale="fit")
  #fig.main %<>% fill_panel(plots1[[2]], row = 2, column=2, scale="fit")
  #fig.main %<>% fill_panel(arrangeGrob(plots1[[3]], plots1[[4]], nrow=2, clip="on"), row = 2, column=3)
  fig.main
  print("Saving Main figure")
  pdf(paste("figures/Fig1-main_", fig, ".pdf", sep=""), width=12.5, height=7.5) # a bit smaller than the figure to avoid padding
    grid.draw(fig.main)
  dev.off()
  
  # supplementary figures show Jaccard+Normalized_counts
  fig.suppl <- grid.arrange(
    arrangeGrob(top=textGrob("A",gp=gpar(fontsize=20,font=2), x=0, hjust=0), 
                plots1[[1]], padding = unit(0, "line")), 
    arrangeGrob(top=textGrob("B",gp=gpar(fontsize=20,font=2), x=0, hjust=0), 
                plots1[[2]], padding = unit(0, "line"), clip="on"), 
    arrangeGrob(top=textGrob("C",gp=gpar(fontsize=20,font=2), x=0, hjust=0),
                plots1[[3]], plots1[[4]], nrow=2, clip="on"),
    arrangeGrob(top=textGrob("D",gp=gpar(fontsize=20,font=2), x=0, hjust=0), 
                plots2[[1]], padding = unit(0, "line")), 
    arrangeGrob(top=textGrob("E",gp=gpar(fontsize=20,font=2), x=0, hjust=0), 
                plots2[[2]], padding = unit(0, "line"), clip="on"), 
    arrangeGrob(top=textGrob("F",gp=gpar(fontsize=20,font=2), x=0, hjust=0),
                plots2[[3]], plots2[[4]], nrow=2, clip="on"),
    layout_matrix=matrix(c(1,2,3,4,5,6), ncol=3, byrow=T))
  
  
  
  #fig.suppl <- multi_panel_figure(width=c(4,4,4), height=c(4,4), unit="inches",
  #                              column_spacing = 0, row_spacing = 0)
  #fig.suppl %<>% fill_panel(plots1[[1]], row = 1, column=1)
  #fig.suppl %<>% fill_panel(plots1[[2]], row = 1, column=2)
  #fig.suppl %<>% fill_panel(arrangeGrob(plots1[[3]], plots1[[4]], nrow=2, clip="on"), row = 1, column=3)
  
  #fig.suppl %<>% fill_panel(plots2[[1]], row = 2, column=1)
  #fig.suppl %<>% fill_panel(plots2[[2]], row = 2, column=2)
  #fig.suppl %<>% fill_panel(arrangeGrob(plots2[[3]], plots2[[4]], nrow=2, clip="on"), row = 2, column=3)
  print("Saving suppl figure")
  pdf(paste("figures/Fig1-suppl_", fig, ".pdf", sep=""), width=12.5, height=8.5) # a bit smaller than the figure to avoid padding
    grid.draw(fig.suppl)
  dev.off()
}

for(fig in names(dd.genes2)){
  multivar_figure(dd.genes2[[fig]], gsub(".RData", "", fig))
}
for(fig in names(dd.exons2)){
  multivar_figure(dd.exons2[[fig]], gsub(".RData", "", fig))
}

#
# OSL/ISL and operon identification
#
find_operons <- function(des){
  # operonic classification
  cl <- rowData(des)$Class
  rowData(des)$Operon <- rep("no TSL", times=length(cl))
  # all ISL genes are downstream genes
  rowData(des)$Operon[cl=="ISL"] <- "downstream"
  # OSL/OSL+ISL genes are either monocistronic or upstream genes
  rowData(des)$Operon[cl=="OSL+ISL" | cl=="OSL"] <- "monocistronic"
  
  # find runs of ISLs on each chromosome and strand
  foreach(chr=iter(unique(rowData(des)$Chr))) %do% {
    for(dd in list(subset(rowData(des), Chr==chr & Strand=="+"), 
                   subset(rowData(des), Chr==chr & Strand=="-"))){
      if(any(dd$Strand=="-")) { dd <- dd[nrow(dd):1,] }
      rl <- rle(dd$Class)
      #for each ISL run, get upstream gene, first ISL and last ISL
      for(rn in which(rl$values=="ISL")){
        isl.end <- sum(rl$lengths[0:rn])
        isl.start <- isl.end-(rl$lengths[rn]-1)
        # upstream gene is right before isl.start
        upstr <- isl.start-1
        
        op.id <- paste("predicted:Operon_", chr, dd$Strand[1], 
                       ifelse(upstr>0, dd$Start[upstr], dd$Start[isl.start]), "-", dd$End[isl.end], sep="")

        print(op.id)
        # label
        #rowData(des)$Operon.annot[rowData(des)$Geneid %in% dd$Geneid[upstr:isl.end]] <- op.id
        #rowData(des)$Operon[rowData(des)$Geneid==dd$Geneid[upstr]] <- "upstream"
        #rowData(des)$Distance[rowData(des)$Geneid==dd$Geneid[isl.end]] <- NA
      }
    }
  }
  return(des)
}

classify_SL <- function(des){
  des <- des[,grep("TSL", colnames(des))]
  cts <- counts(des)
  # combine libraries
  # use geometric mean (enforces at least 1 read per library)

  cts <- sapply(sort(unique(gsub("LIB.*[\\._]", "", colnames(cts)))), function(tsl){
    #rs <- cbind(rowSums(cts[,grep(paste(tsl, "$", sep=""), colnames(cts))]))
    ct <- cts[,grep(paste(tsl, "$", sep=""), colnames(cts))]
    rs <- cbind(apply(ct, 1, function(x) {
        if(any(x)==0) x <- x[-c(which(x==0)[1])]
        exp(mean(log(x)))
      }))
    return(rs)
  })
  rownames(cts) <- rownames(des)
  
  # read counts for OSL and ISL
  OSL <- rowSums(cts[,-c(grep("TSL2|TSL10|TSL12", colnames(cts)))])
  ISL <- rowSums(cts[,grep("TSL2|TSL10|TSL12", colnames(cts))])

  # classify genes by read ratio
  ratio <- OSL/ISL
  cl <- rep("OSL+ISL", times=length(ratio))
  # 0/X = 0. These are ISL genes
  cl[ratio==0] <- "ISL"
  # X/0 = Inf. No 2/10/12 at all. Only OSLs
  cl[is.infinite(ratio)] <- "OSL"
  # 0/0 = NaN. This gene has no TSL coverage. Ignore.
  cl[is.na(ratio)] <- "no TSL"
  # X/X = positive number. The smaller the better (close to 0)
  #cl[ratio<0.2] <- "ISL"
  
  # add to data object
  rowData(des) <- cbind.data.frame(rowData(des), cts)
  rowData(des)$OSL <- OSL
  rowData(des)$ISL <- ISL
  rowData(des)$Ratio <- ratio
  rowData(des)$Class <- cl
  #rowData(des)$Class[! rowData(des)$Expressed] <- "unexpressed"
  return(des)
}

dd.genes2 <- lapply(dd.genes, classify_SL)
dd.exons2 <- lapply(dd.exons, classify_SL)
dd.genes2 <- lapply(dd.genes2, find_operons)
dd.exons2 <- lapply(dd.exons2, find_operons)

#
# write operon predictions as GFF3 files
#
dir.create("operons_GFF3")
write.GFF3 <- function(des, nn=""){
  dd <- subset(rowData(des), !is.na(Operon.annot))
  # cycle through operons to rename them
  operons <- split(dd, dd$Operon.annot)
  
  gff3 <- lapply(1:length(operons), function(i){
    opid <- paste("TSPOP", i, sep="")     # operon name
    gff.1 <- paste(operons[[i]]$Chr[1], ".", "operon", 
              min(operons[[i]]$Start), max(operons[[i]]$End), ".", operons[[i]]$Strand[1], ".", 
              paste("ID=", opid, ";Name=", opid, ";Note=genes:", nrow(operons[[i]]), sep=""), sep="\t")
    gff.2 <- cbind(sapply(1:nrow(operons[[i]]), function(g){
      opgene <- operons[[i]][g,]
      paste(opgene[,"Chr"], ".", "gene", 
            opgene[,"Start"], opgene[,"End"], ".", opgene[,"Strand"], ".", 
            paste("ID=", opid, ".", g, ";Name=", opid, ".", g, " (", gsub(" ", "", opgene[,"Class"]), ");Parent=", opid, sep=""), sep="\t")
    }))
    rbind(gff.1, gff.2)
  })
  gff3 <- rbind("##gff-version 3", paste("# predicted from dataset: ", nn, sep=""), do.call(rbind, gff3))
  write.table(gff3, paste("operons_GFF3/Operons_", nn, ".gff3", sep=""), sep="\t", col.names=F, row.names=F, quote=F)
} 
for(nn in names(dd.genes2)){
  write.GFF3(dd.genes2[[nn]], gsub(".RDS", "", nn))
}
for(nn in names(dd.exons2)){
  write.GFF3(dd.exons2[[nn]], gsub(".RDS", "", nn))
}

#
# make summary tables
#
summ <- data.frame(t(cbind(sapply(dd.genes2, function(x) tapply(colSums(counts(x)), colData(x)$Group2, sum)),
                           sapply(dd.exons2, function(x) tapply(colSums(counts(x)), colData(x)$Group2, sum)))),
                   t(cbind(sapply(dd.genes2, function(x) table(rowData(x)$Class)),
              sapply(dd.exons2, function(x) table(rowData(x)$Class)))),
      Genes=c(sapply(dd.genes2, function(x) nrow(rowData(x))), sapply(dd.exons2, function(x) nrow(rowData(x)))),
      Expressed=c(sapply(dd.genes2, function(x) sum(rowData(x)$Expressed)), sapply(dd.exons2, function(x) sum(rowData(x)$Expressed))),
      Operons=c(sapply(dd.genes2, function(x) length(unique(na.exclude(rowData(x)$Operon.annot)))),
                sapply(dd.exons2, function(x) length(unique(na.exclude(rowData(x)$Operon.annot))))),
      OpGenes=c(sapply(dd.genes2, function(x) sum(!is.na(rowData(x)$Operon.annot))),
        sapply(dd.exons2, function(x) sum(!is.na(rowData(x)$Operon.annot)))))
summ$transpliced <- summ$ISL+summ$OSL+summ$OSL.ISL

Annot <- rep("de novo", times=nrow(summ))
Annot[grep("PRJNA", rownames(summ))] <- "reference"
Stringency <- rep("10 bp", times=nrow(summ))
Stringency[grep("8bp", rownames(summ))] <- "8 bp"
Correction <- rep("raw", times=nrow(summ))
Correction[grep("conv", rownames(summ))] <- "corrected"
cbind(apply(summ, 2, function(x) paste(format(mean(x), big.mark=",", digits=0, scientific = F), format(sd(x), big.mark=",", digits=0, scientific = F), sep=" � ")))
t(aggregate(summ, list(Annot), function(x) paste(format(mean(x), big.mark=",", digits=0, scientific = F), format(sd(x), big.mark=",", digits=0, scientific = F), sep=" � ")))
t(aggregate(summ, list(Stringency), function(x) paste(format(mean(x), big.mark=",", digits=0, scientific = F), format(sd(x), big.mark=",", digits=0, scientific = F), sep=" � ")))
t(aggregate(summ, list(Correction), function(x) paste(format(mean(x), big.mark=",", digits=0, scientific = F), format(sd(x), big.mark=",", digits=0, scientific = F), sep=" � ")))

# statistical tests among groups of datasets
round(apply(summ, 2, function(x) t.test(x~Annot)$p.value), 3)
round(apply(summ, 2, function(x) t.test(x~Stringency)$p.value), 3)
round(apply(summ, 2, function(x) t.test(x~Correction)$p.value), 3)

# operon sizes
lapply(dd.genes2, function(x) table(table(rowData(x)$Operon.annot)))
lapply(dd.exons2, function(x) table(table(rowData(x)$Operon.annot)))
