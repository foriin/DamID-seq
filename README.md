### DamID-seq4NAR
A set of scripts used to analyze DamID-seq experiments in the Laboratory of Analysis of Gene Regulation of IMG RAS, Moscow

Part of these scripts was designed by Ludo Pagie from Bas van Steensel's [research group](https://www.nki.nl/divisions/gene-regulation/van-steensel-b-group/)
from Netherlands Cancer Institute, other part was developed in collaboration with Alexey Pindyurin's [group](https://www.mcb.nsc.ru/laboratory/gatti) from IMCB SB RAS. 
#### Dependencies
To succesfully implement these scripts you need to have these programs installed and running on your system:
  - [Bowtie2](https://github.com/BenLangmead/bowtie2)
  - [Samtools](http://www.htslib.org/)
  - [FASTX-Toolkit](http://hannonlab.cshl.edu/fastx_toolkit/)
  - [Cutadapt](https://github.com/marcelm/cutadapt/tree/77ade52bc2a7fe2d278fdb4256c5b46936011c2c)
  - [HTSeq](https://htseq.readthedocs.io/en/release_0.10.0/)
  
All paths to these programs must be set by yourself in *aligner.sh* and *reads2bins.sh* as bash variables.

R scripts may require packages that are not installed in your library. Please install them from the CRAN repository
or if it doesn't work via [Bioconductor](https://www.bioconductor.org/packages/release/BiocViews.html#___Software).
#### Before you start
If you clone repository you will find make.bins.R script in DamIDseq/bins folder. Please open it, set desirable size
of bins and launch it. It will generate gff and txt files, the former for the HTSeq-count and the latter for
the DamID-seq_analysis.R script. You can re-run this script with different bin size setting and make files for different
bin sizes.

Given that you've set paths to programs and bowtie2 indices in .sh scripts, you need to set path to your fastq files and
path to directory where you want to keep output in parameterfile.txt in DamIDseq folder. (You may change ASSEMBLY parameter if you have indices prepared and know what you're doing). Example:
```
SPECIES=fly
FASTQ_FILES='/home/johndoe/work/projectX/run33/*.gz'
ASSEMBLY=dm3
OUTPUT_DIR='/home/johndoe/work/projectX/OUT'
```

Last but not least, you have to add info about your run files in damid_description.csv file. It's a simple tab-delimited table where in the left column you should write file name of your fastq-file and in the right corresponding info in the format: 
TISSUE.PROTEIN.CONDITIONS.REPLICATE_NUMBER. Example:
```
Data.set  fastq.file
BRAIN.CTCF.vasa(-).1   P155_CGGATG0003.fastq.gz
BRAIN.CTCF.vasa(-).2   P155_ATTGCC0011.fastq.gz
BRAIN.DAM.vasa(-).1   P155_GGCATA0001.fastq.gz
BRAIN.DAM.vasa(-).2   P155_AGTACC0008.fastq.gz
```
