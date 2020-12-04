#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --genome-build ) shift; [[ ! $1 =~ ^- ]] && genome_build=$1;;
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --functional-regions ) shift; [[ ! $1 =~ ^- ]] && functional_regions=$1;;
    --min-repeat-length ) shift; [[ ! $1 =~ ^- ]] && min_repeat_length=$1;;
    --min-repeat-period ) shift; [[ ! $1 =~ ^- ]] && min_repeat_period=$1;;
    --alignments ) shift; [[ ! $1 =~ ^- ]] && alignments=$1;;
    --threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    *) bash ${root}/utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

slop=250
min_coverage=0
max_coverage=200
max_region_length=100000

if [[ "${functional_regions}" == "none" ]]; then
  overlapped_functional_regions=false
else
  overlapped_functional_regions=true
fi

${root}/bin/jq \
  --null-input \
  --arg genome_build ${genome_build} \
  --arg slop ${slop} \
  --arg min_coverage ${min_coverage} \
  --arg max_coverage ${max_coverage} \
  --arg min_repeat_length ${min_repeat_length} \
  --arg min_repeat_period ${min_repeat_period} \
  --arg max_region_length ${max_region_length} \
  --arg functional_regions ${functional_regions} \
  --arg overlapped_functional_regions ${overlapped_functional_regions} \
  '{ 
    "genome build": $genome_build,
    "slop": $slop,
    "minimum coverage": $min_coverage,
    "maximum coverage": $max_coverage,
    "minimum repeat length": $min_repeat_length,
    "minimum repeat period": $min_repeat_period,
    "maximum region length": $max_region_length,
    "overlapped functional regions": $overlapped_functional_regions,
    "functional regions": $functional_regions
  }' \
  > ${output}/make-regions.json

# TODO:
if [[ ! -f ${output}/repeats.bed.gz ]]; then 
  # then download it (given $genome_build),
  # putting the column names in a header line 
  # (see /scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/repeats)
fi 

# TODO
# use (head -1000) to facilitate a small-scale test of tool correctness and tool usage 

filter_repeats_by_length_and_function () {
  bash ${root}/make-regions/filter_by_length.sh \
    --repeats ${repeats} \
    --min-repeat-length ${min_repeat_length} \
    --max-region-length ${max_region_length} \
    --root ${root} \
  | 
  # https://unix.stackexchange.com/a/38311/406037 : 
  if [[ "${functional_regions}" == "none" ]]; then 
    cat
  else 
    bash ${root}/make-regions/filter_by_function.sh \
      --functional-regions ${functional_regions} \
      --root ${root}
  fi 
}

mkdir --parents "${output}/tmp" 
mosdepth_prefix="${output}/tmp/mosdepth.coverage"
${root}/bin/mosdepth \
  --no-per-base \
  --fast-mode \
  --threads ${number_threads} \
  --by <(filter_repeats_by_length_and_function) \
  --fasta ${reference}.fa \
  ${mosdepth_prefix} \
  ${alignments}.cram

${root}/bin/bedtools slop -i ${mosdepth_prefix}.regions.bed.gz -g ${reference}.genome -b ${slop} \
  | awk --assign min_coverage=${min_coverage} '$4 > min_coverage' \
  | awk --assign max_coverage=${max_coverage} '$4 < max_coverage' \
  | awk --assign OFS='\t' '{ print $1, $2, $3 }' \
  | bash ${root}/utilities/sort_compress_index_regions.sh \
    --regions "${output}/regions" \
    --root ${root} \
&& rm -rf "${output}/tmp"

