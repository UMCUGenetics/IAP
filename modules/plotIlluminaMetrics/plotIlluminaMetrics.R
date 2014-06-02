### Plot Illumina picard metrics
### Robert Ernst

# Parse command line arguments
args<-commandArgs(trailingOnly=TRUE)
rootDir = args[1]
runName = args[2]
samples = args[3:length(args)]
#write(samples,file='samples.txt')

### Import library's and plotting functions.
source(paste(rootDir,"plotIlluminaMetrics_include.R",sep="/"))
dir.create("pdfFigures", showWarnings=F) #create output Dir

#Read in hsmetric table
fileName = paste(runName,"HSMetric_summary.txt",sep=".")
hsMetrics = FALSE
colorSet = GetRandomColorSet(length(samples))
if (file.exists(fileName)){
  hsMetrics = TRUE
  summaryTable = read.table(file=fileName, sep="\t", header=TRUE, stringsAsFactors=FALSE)
  summaryTable$sample = as.character(summaryTable$sample)
  summaryTableMelted = melt(summaryTable[,c('sample','PCT_TARGET_BASES_2X','PCT_TARGET_BASES_10X','PCT_TARGET_BASES_20X','PCT_TARGET_BASES_30X','PCT_TARGET_BASES_40X','PCT_TARGET_BASES_50X','PCT_TARGET_BASES_100X')],id.vars = 1)

  #Transpose and write summaryTable
  summaryTableT = t(summaryTable)
  colnames(summaryTableT) = summaryTableT[1,]
  write.table(summaryTableT, file=paste(runName,"HSMetric_summary.transposed.txt", sep="."), col.names=FALSE, na="", quote=FALSE, sep="\t")
  
  pdfOut =  paste("pdfFigures", runName, sep="/")
  
  pctOffBait <- plot_pctOffBait()
  ggsave(paste(pdfOut,"_pctOffBait.pdf", sep=""), pctOffBait, dpi = 300, width=15, height=10)
  
  meanTargetCov <- plot_meanTargetCov()
  ggsave(paste(pdfOut,"_meanTargetCoverage.pdf", sep=""), meanTargetCov, dpi = 300, width=15, height=10)
  
  pctTargetBases <- plot_pctTargetBases()
  ggsave(paste(pdfOut,"_pctTargetBases.pdf", sep=""), pctTargetBases, dpi = 300, width=20, height=10)
}

#Plot metrics to pdf files.
for(i in 1:length(samples)) {
  samplePath = paste(samples[i],"QCStats",samples[i],sep="/")
  
  pdfOut =  paste("pdfFigures", samples[i], sep="/")
  insert_size_metrics = paste(samplePath,"_MultipleMetrics.txt.insert_size_metrics", sep="")
  quality_by_cycle_metrics = paste(samplePath,"_MultipleMetrics.txt.quality_by_cycle_metrics", sep="")
  quality_distribution_metrics = paste(samplePath,"_MultipleMetrics.txt.quality_distribution_metrics", sep="")
  
  insert_size_metrics.table = read.table(file=insert_size_metrics, skip=10, head=TRUE)
  quality_by_cycle_metrics.table = read.table(file=quality_by_cycle_metrics, head=TRUE)
  quality_distribution_metrics.table = read.table(file=quality_distribution_metrics, head=TRUE)
  
  insertSize <- plot_insert_size_metrics()
  ggsave(paste(pdfOut,"_insertSize.pdf", sep=""), insertSize, dpi = 300)
  
  cycleQuality <- plot_quality_by_cycle_metrics()
  ggsave(paste(pdfOut,"_cycleQuality.pdf", sep=""), cycleQuality, dpi = 300)
  
  qualityDistribution <- plot_quality_distribution_metrics()
  ggsave(paste(pdfOut,"_qualityDistribution.pdf", sep=""), qualityDistribution, dpi = 300)
}

### Generate .html based on .Rmd file
workingDir = getwd()
options(knitr.unnamed.chunk.label = runName)
knit(paste(rootDir,"plotIlluminaMetrics.Rmd", sep="/"),quiet=TRUE)
markdownToHTML("plotIlluminaMetrics.md", paste(runName,"picardMetrics.html",sep="."), options=c("use_xhml"), stylesheet=paste(rootDir,"plotIlluminaMetrics_html.css",sep="/"))

### Generate .pdf based on .brew file
brew(paste(rootDir,"plotIlluminaMetrics.brew", sep="/"), "plotIlluminaMetrics.tex")
texi2dvi("plotIlluminaMetrics.tex", pdf = TRUE, clean=TRUE)