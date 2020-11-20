#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    --threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
    --alignments ) shift; [[ ! $1 =~ ^- ]] && alignments=$1;;
    *) bash utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line, 
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }' 

assemble="make-calls/assemble.sh"
regions="${output}/regions"

# get short reads that were originally aligned to given regions

# fetch regions by sweeping through cram (expected to be faster when number of regions > ~10,000) 
# "bedtools intersect -sorted" uses the "chromsweep" algorithm for sorted (-k1,1 -k2,2n) input
# /usr/bin/time --verbose bedtools intersect \
#     -ubam -u -wa -a ${alignments}.cram -b ${regions}.bed -g ${reference}.genome -sorted | 
#   samtools bam2fq > ${regions}.fq

# "random access" of reads (expected to be faster when number of regions < ~10,000)
# "samtools view -M" uses the multi-region iterator 
# (increases speed, removes duplicates and outputs the reads as they are ordered in the file)
/usr/bin/time --verbose samtools view -u -M \
    --threads ${number_threads} \
    -L ${regions}.bed.gz \
    ${alignments}.cram | 
  samtools fastq > ${regions}.fq 

bgzip --force ${regions}.fq 

# generate Makefile for unitig assembly
# https://github.com/lh3/fermikit
# assemble reads from regions into unitigs (-s specifies the genome size and -l the read length)
fermikit_prefix="${output}/fermikit"
make-calls/fermi.kit/fermi2.pl unitig \
    -A ${assemble} \
    -t ${number_threads} \
    -l150 \
    -p ${fermikit_prefix} \
    ${regions}.fq.gz \
  > ${fermikit_prefix}.mak

# execute shell commands to assemble short reads into unitigs
make -f ${fermikit_prefix}.mak || true

[[ ! -e ${output}/filtered_fastq_empty ]] || exit 1

# execute shell commands that align unitigs (with minimap2) and call variants
# ... "-m" means "use minimap2" 
# ... ${fermikit_prefix}.mag.gz contains the (unaligned) unitigs
make-calls/fermi.kit/run-calling \
    -m \
    -t ${number_threads} \
    ${reference}.fa \
    ${fermikit_prefix}.mag.gz | 
  bash -euxo pipefail

