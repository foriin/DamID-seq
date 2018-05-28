#!/usr/bin/R
########################################################################
# Alexey Pindyurin, Anton Ivankin, September 12, 2014
# Artem Ilyin 22.02.18
# For use with bins and no pseudocounts.
# Replicas are just mixed together and go through
# quantile normalization by protein after dam normalization.
# HMM done in this same script
########################################################################

rm(list=ls())

# Load libraries
library(gplots)
library(ggplot2)
library(snapCGH)
library(GenomicRanges)

workDir <- getwd() # or set some different folder
# Declare variables for counting part
###################
prefixDir <- "name_of_your_project" # directory for outputs
binsDir <- file.path(workDir, "/txt_bins")	# location of your GATCs/bins file
sourceDir <- "/path/to/output/dir/from/parameterfile/gatcReadCounts/" # location of HTSeq-count output files. You can specify the highest folder as it is possible. Searching runs recursively.
damIdLocation <- paste0(workDir, "/damid_description.csv") # location of your DamID-Description file
dm.r6 <- F
# Some clarification about format of the damid_description file
# It has to be formatted like this:
#     TISSUE.PROTEIN.conditions.#_of_replica\tname of the fastq.gz file, which was used for mapping and htseq-counting
# e.g. BRAIN.PIWI.vasa(-).1   P155_CGGATG_Piwi1.fastq.gz
#      BRAIN.PIWI.vasa(-).2   P155_ATTGCC_Piwi2.fastq.gz
#      TEST.DAM.wt.1  D1_AGGTTA_LR1_R001.fastq.gz
# and so on

bg4single <- F # do you want to produce bedgraph files for single replicas?
qn <- F # do you want to apply quantile normalization across your tissues?
heatmapColors <- greenred(200)	# color scheme for heatmap

# Declare variables for HMM part

hmm.3 <- F # Set if you want to use three-state HMM instead of two-state
by.chr <- T # Set if you want to apply HMM at chromosomal scale rather than whole genomic
            # scale (default behavior is True)
use.het <- T # Do you want to use heterochromatic loci defined in dm3 Drosophila genome?

################ LOAD FUNCTIONS #################
setwd(workDir)
source("DamID_count.functions.R")
source("DamID_HMM.functions.R")

for (bin in dir(sourceDir))
{
  
  setwd(workDir)
  # Create directories
  foldersNames <- c("Bedgraph", "Statistics", "CSV") # names of folders for different objects created
  outputDirs <- sapply (foldersNames, function(x) file.path(workDir, prefixDir, bin, x))
  lapply(outputDirs, dir.create, showWarnings = FALSE, recursive = T) 
  
  # Make samples list file
  MakeSamplesListFile(file.path(sourceDir, bin), damIdLocation)
  
  
  # Load bins data in data frame
  ################################
  bins.file <- file.path(
    binsDir,
    dir(binsDir)[grepl(paste0(
      "[^0-9]",
      sub("nt", "", bin),
      "[^0-9]"), dir(binsDir))]
  )
  bins <- read.delim(bins.file) %>% arrange(chr, start)
  lastcol <- ncol(bins)
  
  bins <- cbind(bins, matrix(data=NA, nrow=nrow(bins), ncol=nrow(SamplesList)))
  for (i in 1:nrow(SamplesList)){
    colnames(bins)[lastcol+i] <- SamplesList$id[i]
    reads2bins <- fread(SamplesList$path[i], col.names = c("ID", "count")) %>% 
      dplyr::slice(1:(nrow(.) - 5)) %>% mutate(ID = sub("gene:", "", ID),
                                                        chr = sub(".*chr(.*)$", "\\1", ID)) %>% 
      mutate(N = as.integer(sub("^.[a-zA-Z]*(\\d+)", "\\1", chr)),
             chr = sub("^(.[a-zA-Z]*)\\d+", "\\1", chr)) %>% arrange(chr, N)
    if (all(bins$ID == reads2bins$ID)) bins[, lastcol + i] <- reads2bins$count
  }
  rm(i, reads2bins)
  
  
  save.gatc.df <- "Raw.Counts.csv"
  write.table(bins, file=file.path(prefixDir, bin, "CSV", save.gatc.df), sep=";",
              row.names=F, col.names=T, quote=F, dec=".", append=F)
  
  
  
  DATA <- bins
  # Correlation on Counts
  #######################
  MainCorrelations(data=DATA, createPDF=T, corr.desc = "Raw.Counts",
                   file.p = file.path(prefixDir, bin, "Statistics"))
  # Sum replicas
  DATA.s <- DATA[,1:lastcol]
  ids <- unique(gsub("([^.]+\\.)([^.]+\\.)([^.]+\\.).*", "\\1\\2\\3", names(DATA)[(lastcol + 1):ncol(DATA)]))
  for (id in ids){
    DATA.s <- cbind(DATA.s,
                    rowSums(DATA[, grepl(id, names(DATA))]))
    names(DATA.s)[ncol(DATA.s)] <- paste0(id, "sum")
  }
  
  
  
  print("Calculate reads per million")
  
  # Calculation of reads per million
  ###############################
  DATA.rpm <- DATA.s
  
  for (i in (lastcol + 1):(ncol(DATA.rpm))){
      column.sum <- sum(DATA.rpm[, i])
      DATA.rpm[, i] <- DATA.rpm[, i] / column.sum * 10^6
      rm(column.sum)
    }
    rm(i)
    calc.rpm.file <- paste0("RPMs", ".csv")
    write.table(DATA.rpm, file=file.path(prefixDir, bin, "CSV", calc.rpm.file), sep=";",
                row.names=F, col.names=T, quote=F, dec=".", append=F)
  
  # Correlation on Channels
  #########################
  MainCorrelations(data=DATA.rpm, corr.desc = paste0("02.RPMs."),
                   createPDF=T, file.p = file.path(prefixDir, bin, "Statistics"))
  
  # Plot boxplots on RPMs
  ###########################
    bmp(filename=file.path(prefixDir, bin, "Statistics", paste0("03.RPMs.Boxplot", ".bmp")),
        width=2000, height=1000, units="px")
    par(mar=c(12, 8, 0.5, 0.5))
    boxplot(DATA.rpm[, (lastcol + 1):(ncol(DATA.rpm))],
            names=colnames(DATA.rpm)[(lastcol + 1):(ncol(DATA.rpm))], las=2,
            ylab="RPM", outline =F)
    dev.off()
  # Generate BedGraphs if you need to
  if (bg4single) MakeBedGraphFromDATA(data=DATA.rpm,
                                      file.p = file.path(prefixDir, bin, "Bedgraph"))
  
  # DAM Normalization
  ###################
    
  DATA.norm <- DATA.s
  
  DATA.norm <- DATA.norm[, -c((lastcol + 1):ncol(DATA.norm))]
  folks <- names(DATA.s[(lastcol + 1):ncol(DATA.s)])
  uniqueSamples <- unique(gsub("^[^.]+\\.((?:[^.]+\\.)+)sum$", "\\1", folks))
    for (sample in uniqueSamples) {
      folks.sample <- folks[grepl(sample, folks)]
      tissue.id <- folks.sample[!grepl("DAM", folks.sample)]
      dam.id <- folks.sample[grepl("DAM", folks.sample)]
      for (protein in tissue.id) {
        tissue.norm <- paste0(protein, ".norm")
        DATA.norm[[tissue.norm]] <- log2(DATA.rpm[[protein]] / DATA.rpm[[dam.id]])
      }
    }
    for (i in (lastcol + 1):(ncol(DATA.norm))){
      nan.index <- is.nan(DATA.norm[, i])
      inf.index <- is.infinite(DATA.norm[, i])
      DATA.norm[nan.index, i] <- NA
      DATA.norm[inf.index, i] <- NA
      rm(nan.index)
      rm(inf.index)
    }
    rm(i)
  dam.norm <- paste0("Dam.Normalized.csv")
  write.table(DATA.norm, file=file.path(prefixDir, bin, "CSV", dam.norm), sep=";",
              row.names=F, col.names=T, quote=F, dec=".", append=F)
  
  MakeBedGraphFromDATA(data=DATA.norm,
                       file.p = file.path(prefixDir, bin, "Bedgraph"))
  
  # Quantile Normalization
  #########################
  if (qn) {
    data.names <- colnames(DATA.norm)[(lastcol + 1):ncol(DATA.norm)]
    proteins <- unique(gsub("^([^\\.]+)\\..*", "\\1", data.names))
    DATA.o.s.qm <- lapply(proteins, function(prot){
        DATA.norm %>% select(starts_with(prot)) %>% as.matrix() %>% 
          normalize.quantiles() %>% as.data.frame() %>%
          setNames(paste(data.names[grepl(prot, data.names)], "qn", sep = "."))
    })
      
    DATA.norm.qn <- cbind(DATA.norm[, 1:lastcol],
                                     do.call('cbind', DATA.o.s.qm))
    dam.norm.qn <- paste0("Dam.Normalized.QN.",
                          bin, ".csv")
    write.table(DATA.norm.qn, file=file.path(prefixDir, bin, "CSV", dam.norm.qn), sep=";",
                row.names=F, col.names=T, quote=F, dec=".", append=F)
    
    MakeBedGraphFromDATA(data=DATA.norm.qn,
                         file.p = file.path(prefixDir, bin, "Bedgraph"))
  }
  
  
  ####################################################################################
  # HMM
  ####################################################################################
  
  
  
  setwd(file.path(prefixDir, bin))
  dir.create("BioHMM")
  setwd("BioHMM")
  
  # set bin size
  bin.size <- as.integer(sub("nt", "", bin))
  
  # Select chromosomes with which you want to work
  # based on genome release you used and whether
  # you want to add heterochromatin 
  chroms <- c("2L", "2R", "3L", "3R", "X")
  if (!dm.r6 & use.het){
    chroms <- c(chroms, c("2LHet", "2RHet", "3LHet", "3RHet", "XHet", "4", "YHet"))
  }
  
  # Select dataset for HMM based on whether you
  # want to use quantile normalized data or not
  if (qn){
    DATA <- DATA.norm.qn
  } else{
    DATA <- DATA.norm
  }
  
  if (!use.het) DATA <- delete.het(DATA)
  
  if (!by.chr) {
    DATAs <- lapply(DATA[(lastcol + 1):ncol(DATA)], function(x){
    df <- cbind(DATA[1:lastcol], "DamID.value" = x) %>%
      filter(chr %in% chroms, !is.na(DamID.value))
    add <- 0
    for (i in 2:nrow(df)){
      if (df$chr[i - 1] != df$chr[i]) add = df$start[i - 1]
      df$start[i] <- df$start[i] + add
      df$end[i] <- df$start[i] + bin.size
    }
    head(df)
    return(df)
  })
  } else{
    DATAs <- lapply(DATA[(lastcol + 1):ncol(DATA)], function(x){
      names(x) = lapply(chroms, function(y) filter(cbind(DATA[1:lastcol],
                                                         "DamID.value" = x), chr == y,
                                                   !is.na(DamID.value)))
    })
  }
  
  
  if (!hmm.3){
    DATAs.HMM <- hmm.2.wrapper(DATAs, by.chr = by.chr)
  } else {
    DATAs.HMM <- hmm.3.wrapper(DATAs, by.chr = by.chr)
  }
  
  # Prepare data for saving
  DATA.bed <- lapply(DATAs.HMM, function(x){
    bedranges <- x %>% 
      # leave only domains
      filter(domain == 1) %>% 
      # add 'chr' to chromosomes names
      mutate(chrom = paste0("chr", chr))
    # this use of GenomicRanges is to compress resulting bed files and unite extended domains
    gr <- GRanges(seqnames = Rle(bedranges$chrom),  ranges = IRanges(start = bedranges$start, end = bedranges$end)) %>% 
      reduce() 
    # return dataframe which is totally ready to be written in.bed
    data.frame(chr = seqnames(gr), start = start(gr), end = end(gr))
    
  })
  
  
  lapply(seq_along(DATA.bed), function(x){
    tissue = sub("([^.]*)\\.([^.]*)\\..*", "\\2", names(DATA.bed)[x])
    protein = sub("([^.]*)\\.([^.]*)\\..*", "\\1", names(DATA.bed)[x])
    bed.name <- paste0(protein, '.', tissue, '.', bin, '.domains.bed') 
    track_name = paste0('track name="', protein, ' ', tissue, ' HMM"')
    desc = paste0('description="', protein, ' HMM domains for ', tissue, '"')
    write.table(paste(track_name, desc), file=bed.name, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=F)
    write.table(DATA.bed[[x]], file=bed.name, sep="\t", row.names=F, col.names=F, quote=F, dec=".", append=T)
  })
  
}
  
  
  
  
