### Plot Illumina picard metrics
### librarys, plot functions and helper functions.

# Required packages
suppressMessages(library(ggplot2))
suppressMessages(library(knitr))
suppressMessages(library(markdown))
suppressMessages(library(reshape))
suppressMessages(library(xtable))
suppressMessages(library(tools))
suppressMessages(library(brew))

# Plotting functions

#Custom colorscale used for plotting
cbPalette <- c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

plot_insert_size_metrics <- function(){
  ggplot(insert_size_metrics.table, aes(x=insert_size, y=All_Reads.fr_count)) + 
    geom_bar(stat="identity", width=1, fill="#0072B2") +
    xlab("Insert size") + ylab("Count") +
    scale_fill_manual(name="", values=cbPalette)+
    ggtitle(paste("Insert size for all reads in", samples[i], sep=" ")) +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_quality_by_cycle_metrics <- function(){
  ggplot(quality_by_cycle_metrics.table, aes(x=CYCLE, y=MEAN_QUALITY)) + 
    geom_bar(stat="identity", width=1, fill="#0072B2") +
    xlab("Cycle") + ylab("Mean Quality") +
    scale_fill_manual(name="", values=cbPalette)+
    ggtitle(paste("Quality by cycle in", samples[i], sep=" ")) +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_quality_distribution_metrics <- function(){
  ggplot(quality_distribution_metrics.table, aes(x=QUALITY, y=COUNT_OF_Q)) + 
    geom_bar(stat="identity", fill="#0072B2") +
    xlab("Quality Score") + ylab("Observations") +
    scale_fill_manual(name="", values=cbPalette)+
    ggtitle(paste("Quality score distribution in", samples[i], sep=" ")) +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_pctOffBait <- function(){
  ggplot(summaryTable, aes(x=sample, y=PCT_OFF_BAIT, fill=sample)) + 
    geom_bar(stat="identity") +
    xlab("Sample") + ylab("Percentage off bait") +
    scale_fill_manual(name="", values=cbPalette)+
    ggtitle("Percentage off bait") +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          legend.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_meanTargetCov <- function(){
  ggplot(summaryTable, aes(x=sample, y=MEAN_TARGET_COVERAGE, fill=sample)) + 
    geom_bar(stat="identity") +
    xlab("Sample") + ylab("Mean target coverage") +
    scale_fill_manual(name="", values=cbPalette) +
    ggtitle("Mean target coverage") +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          legend.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

plot_pctTargetBases <- function(){
  ggplot(summaryTableMelted,aes(x = sample, y = value)) + 
    geom_bar(aes(fill=variable), stat="identity",position = "dodge") +
    xlab("Sample") + ylab("Percentage") +
    scale_fill_manual(name="", values=cbPalette) +
    ggtitle("Percentage target bases") +
    theme(axis.title = element_text(face="bold", size=15),
          axis.text = element_text(size=15),
          legend.text = element_text(size=15),
          plot.title = element_text(size=15, face ="bold"))
}

# Brew/Tex helper functions
include_graph <- function(width = 1, filename) { 
  paste("\\includegraphics[width=", width, "\\linewidth]{", filename, "}", sep = "") 
}
include_tbl <- function(tableName) {
  print(xtable(tableName), type="latex",table.placement = "H", print.results=FALSE)
}
subfloat_graph <- function(width, filename) { 
  paste("\\subfloat{", "\\begin{minipage}[h]{", width, "\\linewidth}\\centering", include_graph(width = 1, filename), "\\end{minipage}}", sep = "")
}