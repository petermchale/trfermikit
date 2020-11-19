#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --tr_fermikit_output ) shift; [[ ! $1 =~ ^- ]] && tr_fermikit_output=$1;;
    --repeats ) shift; [[ ! $1 =~ ^- ]] && repeats=$1;;
    --functional-regions ) shift; [[ ! $1 =~ ^- ]] && functional_regions=$1;;
    --min-repeat-length ) shift; [[ ! $1 =~ ^- ]] && min_repeat_length=$1;;
    --alignments ) shift; [[ ! $1 =~ ^- ]] && alignments=$1;;
    --number-threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    *) bash utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
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

jq \
  --null-input \
  --arg slop ${slop} \
  --arg min_coverage ${min_coverage} \
  --arg max_coverage ${max_coverage} \
  --arg min_repeat_length ${min_repeat_length} \
  --arg max_region_length ${max_region_length} \
  --arg functional_regions ${functional_regions} \
  '{ 
    slop: $slop,
    min_coverage: $min_coverage,
    max_coverage: $max_coverage,
    min_repeat_length: $min_repeat_length,
    max_region_length: $max_region_length,
    functional_regions: $functional_regions
  }' \
  > ${tr_fermikit_output}/regions.json

filter_repeats_by_length_and_function () {
  bash filter-repeats/filter_by_length.sh \
    --repeats ${repeats} \
    --min-repeat-length ${min_repeat_length} \
    --max-region-length ${max_region_length} \
  | 
  # https://unix.stackexchange.com/a/38311/406037 : 
  if [[ "${functional_regions}" == "none" ]]; then 
    cat
  else 
    bash filter-repeats/filter_by_function.sh ${functional_regions}
  fi 
}

mkdir --parents "${tr_fermikit_output}/tmp" 
mosdepth_prefix="${tr_fermikit_output}/tmp/mosdepth.coverage"
mosdepth \
  --no-per-base \
  --fast-mode \
  --threads ${number_threads} \
  --by <(filter_repeats_by_length_and_function) \
  --fasta ${reference}.fa \
  ${mosdepth_prefix} \
  ${alignments}.cram

bedtools slop -i ${mosdepth_prefix}.regions.bed.gz -g ${reference}.genome -b ${slop} \
  | awk --assign min_coverage=${min_coverage} '$4 > min_coverage' \
  | awk --assign max_coverage=${max_coverage} '$4 < max_coverage' \
  | awk --assign OFS='\t' '{ print $1, $2, $3 }' \
  | bash utilities/sort_compress_index_regions.sh "${tr_fermikit_output}/regions" \
&& rm -rf "${tr_fermikit_output}/tmp"

