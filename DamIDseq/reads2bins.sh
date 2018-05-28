#!/bin/bash
###############################################################################
#
# Ludo Pagie, Tuesday, May 08 2012, reads2GATC_LP120508.sh
# Artem Ilyin, 16.03.18, trying to make it prettier, removing
#	some mostly unnecessary things
#
# DESCRIPTION:
#   simple script to map aligned reads from a bam file onto bins. The
#   bamfile is supplied as an argument to this script.  The genome bins are
#   supplied via CL
#
# ARGUMENTS:
#   - bamfile with mapped reads
#   - gff file with coordinates of bins
#
###############################################################################

if [ $# != 2 ]; then
  echo 'usage: reads2bins.sh mappedReads.bam bins.gff'
  echo 'aborting'
  exit 1
fi

BAMFILE=$1
GFF=$2
# get absolute path of input bam file
SAM_IN=`readlink -e "${BAMFILE}"`

#create links to programs that are going to be used:
SAMT="/path/to/samtools"
HTSEQC="/path/to/htseq-count"
COUNTS_OUT=`basename "${SAM_IN%.bam}_GATCcounts.txt"`

# log some bookkeeping
echo ""
echo "Starting reads2GATC.sh, started as:"
echo "$@"
echo ""
echo "reading from samfile: ${SAM_IN}"
echo ""
echo 'running reads2bins.sh'
echo ''
echo 'using HTSeq-count version:'
echo `${HTSEQC} -h | grep version`
echo ''
echo "using ${GFF} as .gff file"
echo ""


# run HTSeq-count with options:
# -i ID; 'GFF attribute to be used as feature ID'. i
#        The gff file used has features like 'ID=gene:DmelGATCr5chr2L00041'
# -s no; the assay is *not* strand specific
# -o out.sam; name of output file
# -; input from stdin
# ${GFF}; the gff file containing the features to which to map
# -a 0; set minimal alignmentquality to 0, default is 10 for some reason

${SAMT} view -h "${SAM_IN}" | \
  ${HTSEQC} -i ID -s no -a 0 - ${GFF} > ${COUNTS_OUT}
if [ $? -ne 0 ] ; then
  echo "HTSeq-count failed on input file ${SAM_IN}"
  exit 1
fi

exit 0

