### 2. Scripts for mapping and counting reads

Usage: `DamID-seq_count.sh parameterfile.txt bins_size1 bins_size2 ... bins_sizeN`

Please make sure that you have gff files with genomic bins that you're planning to use to count number of reads per bin in the `bins` directory.
If you don't, look into the `make.bins.R` script that is inside `bins` directory.

**Output**: bam files with alignments and read counts from HTSeq-count.
