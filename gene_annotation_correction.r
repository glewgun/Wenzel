# This script corrects gene annotations by splitting genes at
# internal exons that contain distinct peaks of TSL reads.
# 
# The input file is an exon count table obtained via featureCounts.
# The data are processed and stored in DESeq2 format.


library(DESeq2)
library(doParallel)

### 1) load and annotate data from featureCounts
### the need for tidying up chromosome names and other data columns depends 
### on the reference genome and how featureCounts was run

load_fc_data <- function(fl){
  dd <- read.table(fl, header=T, stringsAsFactors = F)
  
  # tidy up sample names
  colnames(dd) <- gsub(".R1.bam", "", gsub("^.*.LIB", "LIB", colnames(dd)))
  
  # tidy up chromosomes
  if(length(grep(";", dd$Chr))>0) dd$Chr <- sapply(strsplit(dd$Chr, ";"), unique)
  
  # tidy up start/end/length
  if(length(grep(";", dd$Start))>0) dd$Start <- sapply(strsplit(dd$Start, ";"), function(x) min(as.numeric(x)))
  if(length(grep(";", dd$End))>0) dd$End <- sapply(strsplit(dd$End, ";"), function(x) max(as.numeric(x)))
  dd$Length <- abs(dd$End-dd$Start)
  
  # tidy up strand info
  dd$Strand[grep("\\+", dd$Strand)] <- "+"
  dd$Strand[grep("\\-", dd$Strand)] <- "-"
  #dd$Strand[grep("\\.", dd$Strand)] <- "."
  print(table(dd$Strand))
  
  # sort by start position (should already be sorted, but to be sure)
  dd <- dd[order(dd$Chr,dd$Start),]
  
  # extract counts matrix
  cts <- as.matrix(dd[,grep("TSL", colnames(dd))])
  rownames(cts) <- dd$Geneid
  
  # if 8bp data, need to lump together TSL 13,14 and 15  
  if(length(grep("8bp", fl))) {
    cts <- cbind(cts[,-grep("TSL13|TSL14|TSL15", colnames(cts))],
                 LIB7135_TSL1345=rowSums(cts[,c("LIB7135_TSL13", "LIB7135_TSL14", "LIB7135_TSL15")]),
                 LIB7136_TSL1345=rowSums(cts[,c("LIB7136_TSL13", "LIB7136_TSL14", "LIB7136_TSL15")]),
                 LIB7137_TSL1345=rowSums(cts[,c("LIB7137_TSL13", "LIB7137_TSL14", "LIB7137_TSL15")]))
  }
  
  # annotate samples with TSL, library and group (TSL2/10/12 or others)
  gr <- rep("others", times=length(colnames(cts)))
  gr[grep("TSL2", colnames(cts))] <- "TSL2"
  gr[grep("TSL10", colnames(cts))] <- "TSL10"
  gr[grep("TSL12", colnames(cts))] <- "TSL12"
  
  gr2 <- gsub("TSL.*", "TSL2/10/12", gr)
  
  ann <- data.frame(ID=colnames(cts), 
                    TSL=gsub("LIB.*_", "", colnames(cts)),
                    LIB=gsub("_.*$", "", colnames(cts)),
                    Group=gr, Group2=gr2)
  # make DESeq2 object
  des <- DESeqDataSetFromMatrix(countData = cts,
                                colData = ann, 
                                design = ~ TSL)
  rowData(des) <- dd[,c(1:6)]
  # how many genes are expressed?
  # use this as "size factor" to adjust percentages later
  gexp <- as.matrix(dd[,grep("LIB", colnames(dd))])
  rowData(des)$Expressed <- rowSums(gexp)>0
  sizeFactors(des) <- sum(rowSums(gexp)>0)
  
  # estimate size factors from background libraries
  #sf <- estimateSizeFactorsForMatrix(as.matrix(dd[,grep(".bam", colnames(dd))]))
  #des$sizeFactor=unname(rep(sf, each=length(des$ID)/3))
  
  print(dim(des))
  return(des)
}

### 2) convert exon-based data to gene-based data
### this is the core function that searches for TSL read peaks at internal exons
### and splits genes at such internal peaks

convert_exons <- function(des){
  cts.thres <- 4	# at least four reads necessary to call a peak
  rowData(des)$Exonid <- rowData(des)$Geneid
  rowData(des)$Geneid <- gsub("exon:", "", gsub("\\.[0-9]*", "", rowData(des)$Geneid))
  rowData(des)$TSLcount <- rowSums(counts(des))
  # only keep counts that occur in all three libraries
  libs <- rowSums(cbind(rowSums(counts(des)[,grep("LIB7135", colnames(counts(des)))]) > 0,
  rowSums(counts(des)[,grep("LIB7136", colnames(counts(des)))]) > 0,
  rowSums(counts(des)[,grep("LIB7137", colnames(counts(des)))]) > 0))
  rowData(des)$TSLcount[libs<3] <- 0
  
  print(length(unique(rowData(des)$Geneid)))

  # convert in parallel across 12 cores
  # process one gene at a time
  cl <- makeCluster(12)
  registerDoParallel(cl)

  print("Converting exons to genes ...")
  splitgenes <- foreach(gene=iter(unique(rowData(des)$Geneid))) %dopar% {
    gene <- subset(des, rowData(des)$Geneid==gene)
    # revert if on negative strand
    if(any(rowData(gene)$Strand=="-")) gene <- gene[nrow(rowData(gene)):1,]
	rowData(gene)$Geneid.new <- rowData(gene)$Geneid
    cts <- rowData(gene)$TSLcount
	# find read peaks along exons
	cts[cts<cts.thres] <- 0
	
	pos <- 0
	if(sum(cts)>0 & length(cts)>1){		# convert only if reads are present
		# find primary peaks
		# defined as an increase in counts from zero at previous exon
		# example exon count vector: 0  0  5  7 12  0  6  0  0 10  4  0  0
		# primary peaks would be: 5, 6 and 10 (vector positions 3, 7 and 10)
		d1 <- diff(c(0,cts))
		pr <- cts>0 & d1-cts==0
		
		# find secondary peaks
		# defined as increase in counts from smaller count at previous exon
		# and decrease in count at following exon.
		# these are an attempt of dealing with difficult (noisy) situations
		# where a single read may span two consecutive exons and was thus assigned to both by featureCounts
		# secondary peaks are only called two exons downstream of primary peaks.
		# This is a conservative way of finding additional peaks and could be improved
		# but would need validation.
		# For example, difficult situations such as 5 7 11 12 can't be resolved at present
		# and are called as a single primary peak at position 1. In contrast, 5 7 12 12 results in
		# a primary peak at position 1 and a secondary peak at position 3. 
		# Such situations should be rare, and the primary peaks alone provide a considerable improvement in any case.
		
		d2 <- c(0,0,pr[0:(length(pr)-2)])	# two-exon shift from primary peaks
		d3 <- c(diff(cts),0)				
		se <- cts>0 & d2==1 & d3<=0
		pos <- which(pr|se)
		
	}
		if(pos>0){  # only check exon positions when peaks are present
		for(g in 1:length(pos)){
			rowData(gene)$Geneid.new[seq(pos[g], length(cts))] <- paste(rowData(gene)$Geneid[1], "-g", g, sep="")
		}  
		# remove exons upstream of the first exon with reads (clearly wrong annotation)
		gene <- gene[grep("-g[0-9]*", rowData(gene)$Geneid.new),]
		}
		# revert back if necessary
		if(any(rowData(gene)$Strand=="-")) gene <- gene[nrow(rowData(gene)):1,]
		# merge exons into genes
		gene.data <- do.call(rbind, lapply(split(rowData(gene), rowData(gene)$Geneid.new), function(x){
		y <- x[1,]
		y$Start <- min(x$Start)
		y$End <- max(x$End)
		y$Length <- abs(y$End-y$Start)
		y$Expressed <- any(x$Expressed)
		return(y)
		}))
		gene.counts <- aggregate(counts(gene), list(rowData(gene)$Geneid.new), sum)
		rownames(gene.counts) <- gene.counts$Group.1
		# for negative strand, the order is wrong way round (factor sorting issue)
		gene.data <- gene.data[match(gene.data$Geneid.new, unique(rowData(gene)$Geneid.new)), ]
		gene.counts <- gene.counts[match(gene.counts$Group.1, unique(rowData(gene)$Geneid.new)), -c(1)]
		return(list(gene.counts, gene.data))
	#} else {
	#	print("No change")
	#	return(list(counts(gene), rowData(gene)))
	#}
  }
  stopCluster(cl)
  print("Merging genes ...")
  # merge all genes into a single DESeq object again
  gc <- do.call(rbind, unname(lapply(splitgenes, function(x){ x[[1]] })))
  gd <- do.call(rbind, unname(lapply(splitgenes, function(x){ x[[2]] })))
  gd$Geneid <- gd$Geneid.new
  if(!identical(rownames(gc), gd$Geneid.new)) print("ERROR during exon conversion!")
  ds <- DESeqDataSetFromMatrix(countData = gc,
                         colData = colData(des),
                         design = ~ TSL)
  rowData(ds) <- gd[, ! colnames(gd) %in% c("Geneid.new", "TSLcount")]
  print(table(sapply(splitgenes, function(x) nrow(x[[2]]))))
  return(ds)
}

### 3) compute distance between genes
get_distances <- function(des){
  for(chr in unique(rowData(des)$Chr)){
  dd <- subset(rowData(des), Chr==chr)
  # + strand
  #print(chr)
  ps <- subset(dd, Strand=="+")  
  di <- c(sapply(1:(length(ps$Start)-1), function(x) { ps$Start[x+1]-ps$End[x] }), NA)
  #print(summary(di))
  #di[di<0] <- 0
  if(nrow(ps)<2) di <- NA
  rowData(des)[rowData(des)$Chr==chr & rowData(des)$Strand=="+", "Distance"] <- di
  # - strand
  ps <- subset(dd, Strand=="-")  
  di <- c(NA, sapply(2:length(ps$Start), function(x) { ps$Start[x]-ps$End[x-1] }))
  #print(summary(di))
  #di[di<0] <- 0
  if(nrow(ps)<2) di <- NA
  rowData(des)[rowData(des)$Chr==chr & rowData(des)$Strand=="-", "Distance"] <- di
  }
  return(des)
}

################### end of functions

args = commandArgs(trailingOnly=TRUE)
fl <- args[1]		# argument is file name of featureCounts count table
# load data
fc <- load_fc_data(fl)

# convert and get distances between genes
fc.conv <- convert_exons(fc)
fc.conv <- get_distances(fc.conv)

# save original and converted datasets as R objects
saveRDS(file=paste(fl, ".RDS", sep=""), fc)
saveRDS(file=paste(fl, ".conv.RDS", sep=""), fc.conv)


