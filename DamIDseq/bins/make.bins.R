# Generate gff file with dm3 genome split onto bins
# for htseq-count and txt file for DamID-seq_analysis.R

# Assuming your working directory is DamIDseq/bins,
# otherwise set it by yourself using setwd()

chroms <- c(chr2L = 23011544,
            chr2LHet = 368872,
            chr2R = 21146708,
            chr2RHet = 3288761,
            chr3L = 24543557,
            chr3LHet = 2555491,
            chr3R = 27905053,
            chr3RHet = 2517507,
            chrX = 22422827,
            chrXHet = 204112,
            chr4 = 1351857,
            chrYHet = 347038)

# SET BIN SIZE
bin.size <- 1000


bins.lst.by.chr <- lapply(names(chroms), function(i){
  bins.st <- seq(1, chroms[[i]], by = bin.size)
  bins.end <- bins.st + bin.size - 1
  data.frame(
    ID = paste0("dm3Bins", bin.size, "nt", i, seq_along(bins.st)),
    chr = i,
    start = bins.st,
    end = bins.end
  )
})

bins.df <- do.call("rbind", bins.lst.by.chr)

write.table(bins.df, paste0("../../Profiles_Domains/txt_bins/bins", bin.size, "nt.dm3.txt"), quote = F, sep = "\t", row.names = F,
            col.names = T)

bins.gff <- data.frame(
  bins.df$chr,
  "src",
  "exon",
  bins.df$start,
  bins.df$end,
  ".",
  ".",
  ".",
  paste0("ID=gene:", bins.df$ID)
)

options(scipen = 999)
write.table(bins.gff, "bins1kdm6.gff", quote = F, sep = "\t", row.names = F,
            col.names = F)
