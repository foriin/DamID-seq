#!/bin/bash
########################################################################
# Ludo Pagie, March 6, 2012, DamID-seq_pipeline_LP120306.sh
# Modified by A. Ilyin for LARG IMG RAS
#
# DESCRIPTION:
#   a pipeline for preprocessing fastq files from DamID-seq experiments
#   the pipeline maps reads onto Drosophila genome dm3.r5 and performs
#   counting by bins of estimated size
#
# DATA:
#   input data is specified in a parameter file which is supplied as
#   argument to this runner script
#
# OUTPUT:
#   data files with read counts per bin in txt and RData format,
#   several tables and plots with QA data
#
########################################################################

usage ()
{
	echo "usage: DamID-seq_count.sh parameterfile.txt bins_size1 bins_size2 ... bins_sizeN"
}

# set code base
CODEDIR=$PWD
ALIGN_SCRIPT="${CODEDIR}/aligner.sh"
READS2GATC_SCRIPT="${CODEDIR}/reads2bins.sh"

ALLSPECIES='human fly'
ALL_ASSEMBLY_HUMAN='hg18 hg19'
ALL_ASSEMBLY_FLY='dm3 dm6'

########################################################################
######## CHECK PARAMETER FILE ##########################################
########################################################################
if [ $# -eq 0 ]
then
echo "no parameterfile given"
usage
exit 1
fi
PARFILE=$1
shift 1 #
if [ ! -f ${PARFILE} ]
then
echo "parfile ${PARFILE} does not exist"
usage
exit 1
else
echo "using par file ${PARFILE}"
fi

# redirect stdout and stderr to a log file
NOW=`date +%d-%m-%Y_%H:%M`
LOG="${PARFILE}_${NOW}.log"
exec &> ${LOG}
# echo some stats to log file
echo "running DamID-seq pipeline"
echo "date = `date`" 
echo "pwd = `pwd`"
echo "commandline ="
echo "$0 $*" 
echo ""

echo "start of initialisation"
echo ""
# import param file
. ./${PARFILE}

# check parameters
# fastq files
if [ -z "${FASTQ_FILES+'xxx'}" ]
then
echo "variable FASTQ_FILES not set, exiting"
exit 1
fi
for fq in ${FASTQ_FILES}
do
if [ ! -e $fq ]
then
echo "file ${fq} does not exist, exiting"
exit 1
fi
done
# SPECIES
if [ -z "${SPECIES+'xxx'}" ]
then
echo "variable SPECIES not set, exiting"
exit 1
fi
correct=0
for sp in ${ALLSPECIES}
do
echo "specs ${sp}"
if [ "${SPECIES}" = "${sp}" ]
then
correct=1
fi
done
if [ ${correct} -eq 0 ]
then
echo "SPECIES should be in ${ALLSPECIES}, exiting"
exit 1
fi
# ASSEMBLY
if [ -z "${ASSEMBLY+'xxx'}" ]
then
echo "variable ASSEMBLY not set, exiting"
exit 1
fi
if [ "${SPECIES}" = 'human' ]
then
ALL_ASSEMBLY=${ALL_ASSEMBLY_HUMAN}
elif [ ${SPECIES} = 'fly' ]
then
ALL_ASSEMBLY=${ALL_ASSEMBLY_FLY}
fi
for as in ${ALL_ASSEMBLY}
do
echo "ass $as"
if [ "${ASSEMBLY}" = "${as}" ]
then
correct=1
fi
done
if [ ${correct} -eq 0 ]
then
echo "for species ${SPECIES} ASSEMBLY should be in ${ALL_ASSEMBLY}, exiting"
exit 1
fi
# OUTPUTDIR
if [ -z "${OUTPUT_DIR+'xxx'}" ]
then
echo "variable OUTPUT_DIR not set, exiting"
exit 1
fi
if [ -d "${OUTPUT_DIR}" ]
then
echo "directory for output exists already: ${OUTPUT_DIR}."
echo "Attempting to continue where previous analysis broke down."
else
mkdir "${OUTPUT_DIR}"
fi

OUTPUT_DIR="`cd \"$OUTPUT_DIR\" 2>/dev/null && pwd || echo \"$OUTPUT_DIR\"`"
echo "end of initialisation"
echo ""

########################################################################
###  END OF CHECK PARAMETER FILE #######################################
########################################################################


#############################
## run bowtie ###############
echo "run bowtie2 for aligning reads to genome"
# check whether we need to restart after failed run
if [ -d ${OUTPUT_DIR}/alignedReads ]; then
	echo "bowtie has been ran succesful previously, skipping this step"
else
	# tempdir for bowtie output
	TMP_DIR="${OUTPUT_DIR}/alignedReads_tmp"
	mkdir ${TMP_DIR}
for df in ${FASTQ_FILES}; do
	D=`dirname "${df}"`
	B=`basename "${df}"`
	A="`cd \"$D\" 2>/dev/null && pwd || echo \"$D\"`/$B"
	ln -s ${A} ${TMP_DIR}
done  
cd ${TMP_DIR}
for fq in ${FASTQ_FILES}; do
	fq_local=`basename ${fq}`
	bash -x ${ALIGN_SCRIPT} ${fq_local} ${ASSEMBLY}

	if [ $? -ne 0 ]; then
		echo "alignscript failed on fastq file; ${fq_local}, aborting" 1>&2
		exit 1
	fi
	rm -f ${fq_local}
done
cd $CODEDIR
# move temp output dir to $OUT_DIR
mv "${TMP_DIR}" "${OUTPUT_DIR}/alignedReads"
fi
echo "end of bowtie2"
echo ""
## end bowtie ###############
#############################

#############################
## run GATC coverage ########
# check whether we need to restart at later stage
echo "run Rscript for read coverage per GATC fragment"
if [ -d ${OUTPUT_DIR}/gatcReadCounts ]; then
	echo "GATC coverage has been run succesful previously, skipping this step"
else
	T_DIR="${OUTPUT_DIR}/gatcReadCounts_tmp"
fi

mkdir ${T_DIR}
mkdir ${OUTPUT_DIR}/gatcReadCounts
set -- "${@/%/nt}" # add nt suffix to all variables that you passed to script starting from 2nd

for fq in ${FASTQ_FILES}; do
	fq_base=`basename ${fq}`
	fq_base=${fq_base%.fastq.gz}
	BAMFILE=${OUTPUT_DIR}/alignedReads/$fq_base/*.bam

	# tempdir for GATC_mapping output
	TMP_DIR="${T_DIR}/$fq_base"

	for nt in "$@" 
	do
		if [ ! -d "${TMP_DIR}" ]; then
			mkdir ${TMP_DIR}
		fi
		cd ${TMP_DIR}
		GFF=`find "${CODEDIR}/bins/" -regex ".*[^0-9]${nt%nt}[^0-9].*"`
		echo "${nt} ${GFF}"
		if [ ! -f "${GFF}" ]; then
			echo "No gff of this bins size found"
			continue
		fi
		bash -x ${READS2GATC_SCRIPT} ${BAMFILE} ${GFF}
		if [ $? -ne 0 ]; then
			echo "bins coverage script failed on ${BAMFILE}, aborting"
			exit 1
		fi
		cd ${CODEDIR}
		if [ ! -d "${OUTPUT_DIR}/gatcReadCounts/${nt}" ]; then
			mkdir ${OUTPUT_DIR}/gatcReadCounts/${nt}
		fi
		mv -f ${TMP_DIR} ${OUTPUT_DIR}/gatcReadCounts/${nt}/$fq_base
		if [ $? -ne 0 ]; then
			echo "cannot move temp output dir to ${OUTPUT_DIR}/gatcReadCounts/${nt}/$fq_base, aborting" 1>&2
			exit 1
		fi
	done
	
done

rm -R ${T_DIR}
