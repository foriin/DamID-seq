####################################
#  Make samples list file function #
####################################


MakeSamplesListFile <- function(SOURCE, DAMID){
  filePath <- list.files(path=SOURCE, pattern=paste0("*", "_local_GATCcounts.txt"), full.names=T, recursive=T)
  # baseFile <- unique(sub("(.*)_(edge|inner).*", "\\1", basename(filePath), perl=T))
  damid.info <- read.delim(DAMID, header=T, sep="\t", stringsAsFactors=F)
  info.split <- sapply(filePath, function(fp){
    print(sub("([0-9_.a-zA-Z-]+)_local_GATCcounts.txt", "\\1",
              basename(fp)))
    descriptor <- damid.info[,1][grepl(sub("([0-9_.a-zA-Z-]+)_local_GATCcounts.txt", "\\1",
                                           basename(fp)), damid.info[,2])]
    print(descriptor)
    subs <- c(paste0("\\2.\\1.\\3.\\4"), "\\2", "\\1", "\\3", "\\4")
    
    sapply(subs, function(x) cbind(sub("([0-9_a-zA-Z-]+)\\.([0-9_a-zA-Z-]+)\\.([0-9_a-zA-Z-]+)\\.([0-9]+)", x, descriptor)))
  })
  SamplesList <<- data.frame(unname(t(info.split)), filePath, stringsAsFactors = F)
  names(SamplesList) <<- c("id", "protein", "tissue", "cond", "rep", "path")
}

######################################
# Functions for statistics retrieval #
######################################

MakeCorMat <- function(data, cormethod, use.vals = "complete.obs", lastcol = 4){
  cormat <- as.matrix(data[,c((lastcol+1):ncol(data))]) %>%
    cor(method = cormethod, use = use.vals) %>% 
    round(digits = 2)
}


# Main correlations function
############################
MainCorrelations <- function(data, use.opt="everything", corr.desc, createPDF=T, file.p, ...) {  
  lapply(c("spearman", "pearson"), function(meth){
    cors <- MakeCorMat(data, meth)
    assign(paste0(corr.desc, ".", meth, ".cor"), cors)
    if (createPDF == T){
      options(warn=-1)
      pdf(file=file.path(file.p, paste0(corr.desc,  "_", meth, "_correlation_heatmap", ".pdf")), width=12, height=12)
      heatmap.2(cors, col=heatmapColors, breaks=seq(from=-1, to=1, by=0.01),
                Rowv=T, Colv=T, dendrogram="both", trace="none", cellnote=cors,
                notecol="white", notecex = 0.5, margins=c(7,7),
                main=paste(meth, "'s correlation coefficients and hierarchical clustering for\n", "'",
                           sub("\\d+\\.(.*)", "\\1", corr.desc), "'", sep=""),
                cex.lab=1.1, cexRow=0.6, cexCol=0.6, lmat=matrix(c(4,2,3,1), ncol=2),
                lwid=c(0.1, 0.9), lhei=c(0.15, 0.85), key=T, density.info="density")
      options(warn=0)
      dev.off()
    }  
  })
}

# ACF on data function
######################
AcfOnData <- function(data, method, acf.desc, lastcol = 4, ylab.val, na.data, file.p) {
  pdf(file=file.path(file.p, paste(acf.desc, "_", "ACF_plot.pdf", sep="")), width=11.69, height=8.27)
  par(mfrow=c(3, 4))
  par(mai=c(0.7, 0.7, 0.7, 0.5))
  DATAs.acf <- data
  descr <- "all"
  acf.order.list <- order(names(DATAs.acf)[(lastcol + 1):ncol(DATAs.acf)]) + lastcol
  if (method == "acf") {    
    for (i in acf.order.list) {
      if (na.data == T) {
        acf.na.data <- sum(!is.na(DATAs.acf[, i]))
      } else {
        acf.na.data <- sum(DATAs.acf[, i] != 0)
      }
      acf(DATAs.acf[, i], na.action=na.pass, main=paste(names(DATAs.acf[i]), "\n(", descr," GATCs: ", acf.na.data, "\nout of ", nrow(DATAs.acf), ")", sep=""), ylab=ylab.val)
    }
  } else {
    for (i in acf.order.list) {
      if (na.data == T) {
        acf.na.data <- sum(!is.na(DATAs.acf[, i]))
      } else {
        acf.na.data <- sum(DATAs.acf[, i] != 0)
      }
      plot(density(DATAs.acf[, i], na.rm=T), main=paste(names(DATAs.acf[i]), "\n(", descr," GATCs: ", acf.na.data, "\nout of ", nrow(DATAs.acf), ")", sep=""), ylab=ylab.val)
    }
  }		
  rm(i)
  dev.off()
  
}


ScatCor <- function(data, lastcol = 4, pref, file.p){
  png.name <- paste0(pref, ".Scatter_Plots_and_Correlations.png")
  corplot <- ggpairs(
    data[, (lastcol + 1):ncol(data)],
    title = "Scatter Plots and Pearson Correlations",
    upper = list(
      continuous = wrap("cor", size = 15)),
    lower = list(
      continuous=wrap("smooth", colour="blue")
    ),
    diag = NULL) +
    theme_grey(base_size = 20)
  png(filename = file.path(file.p, png.name), width = 1000, height = 1000)
    print(corplot)
  dev.off()
}


# Function for calculating mean values of logarithms
log2Mean <- function(x, y){
  log2((2^x + 2^y)*0.5)
}


############################
# VISUALIZATION/PROFILES   #
############################

MakeBedGraphFromDATA <- function(data, lastcol = 4, sel.cols = 2:4, file.p) {
  options(scipen = 999)
  data.bg <- data
  #add chr to the names of chromosomes
  data.bg$chr <- paste("chr", data.bg$chr, sep='')
  # Calculate Bedgraph
  for (j in (lastcol + 1):(ncol(data.bg))) {
    bg.file <- file.path(file.p, paste0(names(data.bg)[j], ".bedgraph"))
    selected.data <- data.bg[, c(sel.cols, j)]
    # selected.data[, 3] <- selected.data[, 3] - 3
    selected.data <- selected.data[!is.na(selected.data[, ncol(selected.data)]), ]
    write.table(paste0("track type=bedGraph name='", names(data.bg)[j], ".",
                       "' description='bedgraph_profile' visibility=full autoScale=on alwaysZero=on"),
                file=bg.file, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
    write.table(selected.data, file=bg.file, sep=" ", row.names=F,
                col.names=F, quote=F, dec=".", append=T)
  }
  rm(j)		
  
}