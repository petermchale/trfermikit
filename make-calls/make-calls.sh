#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    *) bash ${root}/utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
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

read_config () { 
  local key1_=$1 
  local key2_=$2
  ${root}/utilities/read_config.sh ${root} ${output} ${key1_} ${key2_}
}

reference=$(read_config general reference)
number_threads=$(read_config general numberThreads)
alignments=$(read_config general alignments) 

single_base_match_reward=$(read_config makeCalls singleBaseMatchReward)
single_base_mismatch_penalty=$(read_config makeCalls singleBaseMismatchPenalty)
gap_open_penalties=$(read_config makeCalls gapOpenPenalties)
gap_extension_penalties=$(read_config makeCalls gapExtensionPenalties)
minimum_unitig_mapping_quality=$(read_config makeCalls minUnitigMappingQuality)

# we will get short reads that were originally aligned to the following regions
regions="${output}/regions"

# fetch regions by sweeping through cram (expected to be faster when number of regions > ~10,000) 
# "bedtools intersect -sorted" uses the "chromsweep" algorithm for sorted (-k1,1 -k2,2n) input
# /usr/bin/time --verbose bedtools intersect \
#     -ubam -u -wa -a ${alignments}.cram -b ${regions}.bed -g ${reference}.genome -sorted | 
#   bin/samtools bam2fq > ${regions}.fq

# "random access" of reads (expected to be faster when number of regions < ~10,000)
# "samtools view -M" uses the multi-region iterator 
# (increases speed, removes duplicates and outputs the reads as they are ordered in the file)
/usr/bin/time --verbose ${root}/bin/samtools view -u -M \
    --threads ${number_threads} \
    -L ${regions}.bed.gz \
    ${alignments}.cram | 
  ${root}/bin/samtools fastq |
  ${root}/bin/bgzip --stdout \
  > ${regions}.fq.gz 
# NOTE: the above code pulls down unmapped reads with RNAME and POS contained in ${regions}.bed.gz,
# as can be checked by running "samtools view -h ${alignments}.cram ${region} | samtools view -f 4"

# generate Makefile for unitig assembly
# https://github.com/lh3/fermikit
# assemble reads from regions into unitigs (-s specifies the genome size and -l the read length)
do_assembly="${root}/make-calls/assemble.sh"
fermikit_prefix="${output}/fermikit"
assembly_diagnostics="${fermikit_prefix}.assembly.diagnostics.json"
${root}/make-calls/fermi.kit/fermi2.pl unitig \
    -A ${do_assembly} \
    -d ${assembly_diagnostics} \
    -x ${root} \
    -t ${number_threads} \
    -l150 \
    -p ${fermikit_prefix} \
    ${regions}.fq.gz \
  > ${fermikit_prefix}.mak

# execute shell commands to assemble short reads into unitigs
make -f ${fermikit_prefix}.mak || true

filtered_fastq_empty () {
  # echo $(${root}/bin/jq --raw-output '.filteredFastqEmpty' ${assembly_diagnostics})
  ${root}/bin/jq --raw-output '.filteredFastqEmpty' ${assembly_diagnostics}
}

if [[ $(filtered_fastq_empty) == "true" ]]; then 
  exit 1
fi

# execute shell commands that align unitigs (with minimap2) and call variants
# ... "-m" means "use minimap2" 
# ... ${fermikit_prefix}.mag.gz contains the (unaligned) unitigs
${root}/make-calls/fermi.kit/run-calling \
    -m \
    -t ${number_threads} \
    -A ${single_base_match_reward} \
    -B ${single_base_mismatch_penalty} \
    -O ${gap_open_penalties} \
    -E ${gap_extension_penalties} \
    -q ${minimum_unitig_mapping_quality} \
    -r ${root} \
    ${reference}.fa \
    ${fermikit_prefix}.mag.gz \
  | bash -euxo pipefail

calls="${output}/fermikit.raw"
gunzip --force ${calls}.vcf.gz 
cat ${calls}.vcf \
  | bash ${root}/utilities/sort_compress_index_calls.sh \
    --calls ${calls} \
    --root ${root}


