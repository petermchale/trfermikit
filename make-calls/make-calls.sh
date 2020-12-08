#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    --threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
    --alignments ) shift; [[ ! $1 =~ ^- ]] && alignments=$1;;
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

single_base_match_reward="10" 
single_base_mismatch_penalty="12" 
gap_open_penalties="6,26" # there are two because the cost function of gap length is piecewise linear
gap_extension_penalties="1,0" # there are two because the cost function of gap length is piecewise linear
minimum_unitig_mapping_quality="1" 

${root}/bin/jq \
  --null-input \
  --arg single_base_match_reward ${single_base_match_reward} \
  --arg single_base_mismatch_penalty ${single_base_mismatch_penalty} \
  --arg gap_open_penalties ${gap_open_penalties} \
  --arg gap_extension_penalties ${gap_extension_penalties} \
  --arg minimum_unitig_mapping_quality ${minimum_unitig_mapping_quality} \
  '{ 
    "single-base match reward": $single_base_match_reward,
    "single-base mismatch penalty": $single_base_mismatch_penalty,
    "gap-open penalties": $gap_open_penalties,
    "gap-extension penalties": $gap_extension_penalties,
    "minimum unitig mapping quality": $minimum_unitig_mapping_quality
  }' \
  > ${output}/make-calls.json

# get short reads that were originally aligned to given regions
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
  ${root}/bin/samtools fastq > ${regions}.fq 

${root}/bin/bgzip --force ${regions}.fq 

# generate Makefile for unitig assembly
# https://github.com/lh3/fermikit
# assemble reads from regions into unitigs (-s specifies the genome size and -l the read length)
do_assembly="${root}/make-calls/assemble.sh"
fermikit_prefix="${output}/fermikit"
assembly_diagnostics="${fermikit_prefix}.assembly.diagnostics.json"
${root}/make-calls/fermi.kit/fermi2.pl unitig \
    -A ${do_assembly} \
    -d ${assembly_diagnostics} \
    -r ${root} \
    -t ${number_threads} \
    -l150 \
    -p ${fermikit_prefix} \
    ${regions}.fq.gz \
  > ${fermikit_prefix}.mak

# execute shell commands to assemble short reads into unitigs
make -f ${fermikit_prefix}.mak || true

filtered_fastq_empty () {
  echo $(${root}/bin/jq --raw-output '."filtered fastq empty"' ${assembly_diagnostics})
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


