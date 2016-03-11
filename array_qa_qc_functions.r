library(gdata)
library(oligo)
library(pd.mogene.2.1.st)
library(mogene21sttranscriptcluster.db)
library(Heatplus)
library(ggplot2)
library(reshape2)

## Help Documentation
describe = function(obj) {
  if ('help' %in% names(attributes(obj))) {
    writeLines(attr(obj, 'help'))
  }
}
attr(describe, 'help') = "
This function prints the contents of the 'help' attribute of any R object. 
It is meant to provide help documentation in the same vein as Docstrings in Python. 
"

make.colors.boxplot <- function(col.strains)
{
    use.cols <- darkColors(length(unique(col.strains)))
    names(use.cols) <- unique(col.strains)
    return(as.character(use.cols[col.strains]))
}

make.boxplot <- function(use.exprs, type=c("raw", "norm"), make.pdf=F, base.name="test", order.by=c("Mating", "Sex", "Number"), color.by="Mating", highlight.names=NULL, ...)
{
  type <- match.arg(type)
  
  #use.exprs <- switch(type, raw=rma(raw.exprs, normalize=FALSE, target="core"), norm=rma(raw.exprs, normalize=TRUE, target="core"))
  
  basic.ord <- do.call("order", pData(use.exprs)[,order.by,drop=F])
  
  use.exprs <- use.exprs[,basic.ord]
  
  layout(c(1,1,1,1,2,3,3,3))
  par(mar=c(0,4,4,2))
  
  if (type == 'raw') {title = "Raw Expression (Background Corrected)"}
  if (type == 'norm') {title = "Normalized Expression"}
  boxplot(exprs(use.exprs), xaxt="n", ylab="Expression", main=title, outline=FALSE, col=make.colors.boxplot(pData(use.exprs)[,color.by]))
  
  #RIN subplot
  if ('RIN' %in% names(pData(use.exprs)))
  {
    par(mar=c(0,4,1,2))
    pData(use.exprs)$RIN = as.numeric(as.character(pData(use.exprs)$RIN))
    plot(x=1:ncol(use.exprs), y=pData(use.exprs)$RIN, xlim=c(1, ncol(use.exprs)) + c(-.5, .5), ylim=range(pData(use.exprs)$RIN, na.rm=T), xaxt="n", type="p", ylab="RIN", xlab="", col=make.colors.boxplot(pData(use.exprs)[,color.by]))
    abline(h=7, lty="dashed")
  } else {
	warning("RIN values were not found in the sample annotations.")
  }
  
  #labels and such
  if (missing(highlight.names) || is.null(highlight.names))
  {
    Axis(at=1:ncol(use.exprs), side=1, labels=colnames(use.exprs), las=2, ...)
  }
  else if(is.character(highlight.names) && all(highlight.names %in% colnames(use.exprs)))
  {
    should.col <- ifelse(colnames(use.exprs) %in% highlight.names, 'red', 'black')
    rle.col <- Rle(should.col)
    
    lo.pos <- 1
    
    for (i in 1:nrun(rle.col))
    {
      positions <- lo.pos:(runLength(rle.col)[i] + lo.pos-1)
      Axis(at=positions, side=1, labels=colnames(use.exprs)[positions], las=2, col.axis=runValue(rle.col)[i], ...)
      lo.pos <- max(positions) + 1
    }
    
    Axis(at=1:ncol(use.exprs), side=1, labels=colnames(use.exprs), las=2, col.axis=ifelse(should.col, 'red', 'black'), ...)
    
  }
  
  if (make.pdf) {
    d = dev.copy2pdf(file=paste0(base.name, '_', type, '_boxplot.pdf'), width=16, height=10)
	writeLines(paste0("Figure saved to: ", paste0(base.name, '_', type, '_boxplot.pdf')))
	d = dev.off()
  }
}
attr(make.boxplot, 'help') = "
This function creates a boxplot of the raw or normalized expression values. If RIN values 
are available a subplot will be added.

Parameters:
use.exprs: An ExpressionSet returned by the rma() function.
type: The type of plot to create, either 'raw' or 'norm'. 
make.pdf: A logical indicating if a PDF should be created. 
base.name: The base filename of the PDF to create (default is 'test') 
order.by: A character vector containing the column names that will be used to order 
   the samples (default is c('Mating', 'Sex', 'Number')).
color.by: The column name used to assign colors to the samples (default is 'Mating'). 
highlight.names: A character vector containing sample names (default=NULL). Can be used 
   to highlight samples in the plot. 
...: Additional parameters can be passed to format the x-axis labels.
"


#A convienience function for making bacterial spike plots
#Input:
#       raw.exprs: A GeneFeatureSet or similar object which has an rma method returning an ExpressionSet
#       base.name: The base name of the PDF to generate (if make.pdf == T)
#       pgf.file: probe group file from Affy for the appropriate array type.
#Output:
#       A PDF named as base.name_bac_spikes.pdf

plot.bac.spikes <- function(use.exprs, pgf.file, base.name="test", make.pdf=F)
{
    require(ggplot2)
    require(affxparser)
    #use.exprs <- rma(raw.exprs)
    
	if (file.exists(pgf.file)) {
		pgf.list <- readPgf(pgf.file)
	} else {
		stop("Probe group file not found!")
	}
    probeset.types <- data.frame(fsetid=pgf.list$probesetId, fsetName=pgf.list$probesetName, type=pgf.list$probesetType, stringsAsFactors=FALSE)
    
    spike.probes <- probeset.types[probeset.types$type == 'control->affx->bac_spike',]

    spike.exprs <- exprs(use.exprs)[as.character(unique(spike.probes$fsetid)),]
    
    spike.dta <- melt(spike.exprs)
    spike.dta.merge <- merge(spike.dta, spike.probes, by.x="Var1", by.y="fsetid", all=F, incomparables=NA, sort=FALSE)
    
    spike.dta.merge$`Bac. Spike` <- sapply(strsplit(spike.dta.merge$fsetName, "-"), "[[", 4)
    #should be BioB<BioC<BioD<Cre
    spike.dta.merge$`Bac. Spike` <- factor(spike.dta.merge$`Bac. Spike`, levels=c("cre", "bioD", "bioC", "bioB"), labels=c("Cre", "BioD", "BioC", "BioB"), ordered=TRUE)
    
    bac.plot <- qplot(x=Var2, y=value, group=`Bac. Spike`, color=`Bac. Spike`, data=spike.dta.merge, stat="summary", fun.y=mean, geom="line", ylab="log2(Expression)", xlab="", main="Bacterial Spikes")
	bac.plot <- bac.plot + theme(axis.text.x=element_text(size=8, angle=90, hjust=1))
    
	plot(bac.plot)
	
    if (make.pdf) {
        dev.copy2pdf(file=paste0(base.name, "_bac_spikes.pdf"), width=16, height=10)
		writeLines(paste0("Figure saved to: ", paste0(base.name, "_bac_spikes.pdf")))
		d = dev.off()
    }
}
attr(plot.bac.spikes, 'help') = "
This function creates a bacterial spike plot.

Parameters:
use.exprs: An ExpressionSet returned by the rma() function.
pgf.file: A probe group file for the array.
make.pdf: A logical indicating if a PDF should be created. 
base.name: The base filename of the PDF to create (default is 'test').
"

#A convienience function for making polya spike plots
#Input:
#       raw.exprs: A GeneFeatureSet or similar object which has an rma method returning an ExpressionSet
#       base.name: The base name of the PDF to generate (if make.pdf == T)
#       plot.type: Either "AFFX-r2-Bs" or "AFFX" which are the two types of probesets, defaults to "AFFX-r2-Bs"
#       pgf.file: probe group file from Affy for the appropriate array type.
#Output:
#       A PDF named as base.name_polya_spikes.pdf
plot.polya.spikes <- function(use.exprs, pgf.file, plot.type=c("AFFX-r2-Bs", "AFFX"), base.name="test", make.pdf=F)
{
  require(ggplot2)
  require(affxparser)
  
  plot.type <- match.arg(plot.type)
  
  #use.exprs <- rma(raw.exprs)
  
  pgf.list <- readPgf(pgf.file)
  probeset.types <- data.frame(fsetid=pgf.list$probesetId, fsetName=pgf.list$probesetName, type=pgf.list$probesetType, stringsAsFactors=FALSE)
  
  poly.a.spikes <- probeset.types[probeset.types$type == "control->affx->polya_spike",]
  
  spike.exprs <- exprs(use.exprs)[as.character(unique(poly.a.spikes$fsetid)),]
  
  spike.dta <- melt(spike.exprs)
  spike.dta.merge <- merge(spike.dta, poly.a.spikes, by.x="Var1", by.y="fsetid", all=F, incomparables=NA, sort=FALSE)
  
  spike.dta.merge$type <- ifelse(grepl("AFFX-r2-Bs", spike.dta.merge$fsetName), "AFFX-r2-Bs", "AFFX")
  spike.dta.merge$pos <- sub("_*s*_[as]t", "", sapply(strsplit(spike.dta.merge$fsetName, "-"), function(x) x[length(x)]))
  spike.dta.merge$pos <- factor(spike.dta.merge$pos, levels=c("5", "M", "3"), ordered=T)
  
  spike.dta.merge$direction <- sapply(strsplit(spike.dta.merge$fsetName, "[-_]"), function(x) x[length(x)])
  
  
  split.fset <- strsplit(spike.dta.merge$fsetName, "-")
  spike.dta.merge$Spike <- ifelse(spike.dta.merge$type == "AFFX", sapply(split.fset, "[", 2), sapply(split.fset, "[", 4))
  spike.dta.merge$Spike <- sub("X", "", spike.dta.merge$Spike)
  substr(spike.dta.merge$Spike, 1, 1) <- toupper(substr(spike.dta.merge$Spike, 1, 1))
  
  #from the affy exon/gene array whitepaper: lys<phe<thr<dap, not sure about Trpn
  spike.dta.merge$Spike <- factor(spike.dta.merge$Spike, levels=c("Dap", "Thr", "Phe", "Lys", "Trpn"), ordered=T)
  
  pa.basic.plot <- qplot(x=Var2, y=value, group=Spike, color=Spike, data=spike.dta.merge[spike.dta.merge$type == plot.type,], stat="summary",
                         fun.y=mean, geom="line", ylab="log2(Expression)", xlab="", main="PolyA Spikes", facets=pos~direction) + theme(axis.text.x=element_text(size=8, angle=90, hjust=1))
  
  
  plot(pa.basic.plot)
	
  if (make.pdf) {
	dev.copy2pdf(file=paste0(base.name, "_polya_spikes.pdf"), width=16, height=10)
	writeLines(paste0("Figure saved to: ", paste0(base.name, "_polya_spikes.pdf")))
	d = dev.off()
  }
}
attr(plot.polya.spikes, 'help') = "
This function creates a polyA spike plot.

Parameters:
use.exprs: An ExpressionSet returned by the rma() function.
pgf.file: A probe group file for the array.
plot.type: A string indicating the type of control probesets for the array (default is 'AFFX-r2-Bs').
make.pdf: A logical indicating if a PDF should be created. 
base.name: The base filename of the PDF to create (default is 'test').
"


#A function that provides a consistent way of processing expression data for heatmaps etc.
#Input:
#   raw.exprs needs to be FeatureSet
#   num.genes should either be an integer value indicating the 'num.genes' most variable genes to return
#           or should not be specified at all.       
#Output:
#   An ExpressionSet
heatmap.process <- function(use.exprs, num.genes)
{
    norm.exprs.basic <- use.exprs
    norm.exprs.basic.mat <- exprs(norm.exprs.basic)
    if(missing(num.genes) || is.null(num.genes) || is.na(num.genes))
    {
        return(norm.exprs.basic.mat)
    } else {
		exprs.var <- apply(norm.exprs.basic.mat, 1, var)
        norm.exprs.basic.var <- norm.exprs.basic.mat[order(exprs.var, decreasing=TRUE),][1:num.genes,]
        return(norm.exprs.basic.var)
    }
   
}
attr(heatmap.process, 'help') = "
This function processes expression data for heatmaps.

Parameters:
use.exprs: An ExpressionSet returned by the rma() function.
num.genes: A number indicating the number of most variable genes to include. If NULL, all genes
   will be included (default is 1000).

Returns:
A matrix of expression values.
"

make.heatmap <- function(use.exprs, cut.dist=NULL, num.genes=1000, base.factors=c("Sex"), rin.breaks=c(0, 5, 7, 10), make.pdf=F, base.name="test")
{
  
  norm.exprs.basic.var <- heatmap.process(use.exprs, num.genes)
  
  all.dta <- pData(use.exprs)
  
  stopifnot(all(rownames(all.dta) == colnames(norm.exprs.basic.var)))
  
  if (missing(cut.dist) || is.null(cut.dist) || is.na(cut.dist))
  {
    cluster <- NULL
  } else {
    cluster <- list(cuth=cut.dist)
  }
  
  ## Add option of plotting RIN as continuous variable
  rin.discrete=T
  if (rin.discrete) {
	if (!is.null(rin.breaks) && 'RIN' %in% names(all.dta)) {
		all.dta$RIN = as.numeric(as.character(all.dta$RIN))
		all.dta$RIN <- cut(all.dta$RIN, rin.breaks)
		use.dta <- all.dta[,append(base.factors, "RIN", after=0)]
	} else {
		warning("RIN values were not found in the sample annotations.")
		use.dta <- all.dta[,base.factors]
	}
  } else {
    all.dta$RIN = as.numeric(as.character(all.dta$RIN))
	use.dta <- all.dta[,append(base.factors, "RIN", after=0)]
  }
  
  new.heat <- annHeatmap(norm.exprs.basic.var, annotation=use.dta, dendrogram = list(clustfun = hclust, distfun = dist, Col = list(status = "yes"), Row = list(status = "hidden")),
                         labels=list(Row=list(labels=rep("", nrow(norm.exprs.basic.var))), Col=list(nrow=0, labels=rep("", ncol(norm.exprs.basic.var)))), legend = TRUE, cluster = cluster)
  
  # write colnames of clustered IDs to file
  #write.table(file="./bat_clustered.txt",x=colnames(new.heat$data$x2),quote=F)
  
  #cuth was chosen by examination of plot(new.heat$dendrogram$Col$dendro)
  
  plot(new.heat)
  
  if (make.pdf) {
      d = dev.copy2pdf(file=paste0(base.name, '_heatmap.pdf'), width=16, height=10)
      writeLines(paste0("Figure saved to: ", paste0(base.name, '_heatmap.pdf')))
	  d = dev.off()
  }
}
attr(make.heatmap, 'help') = "
This function creates an annotated heatmap for expression data.

Parameters:
use.exprs: An ExpressionSet returned by the rma() function.
cut.dist: The height at which to cut the dendrogram (default is NULL; no cutting).
num.genes: A number indicating the number of most variable genes to include. If NULL, all genes
   will be included (default is 1000).
base.factors: A character vector containing column names that will be used to annotate the
   heatmap (default is c('Sex')).
rin.breaks: A numeric vector indicating the break points for discretizing the RIN 
   scores (default is c(0, 5, 7, 10)).
make.pdf: A logical indicating if a PDF should be created. 
base.name: The base filename of the PDF to create (default is 'test').
"