#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --repeats ) shift; [[ ! $1 =~ ^- ]] && repeats=$1;;
    --min-repeat-length ) shift; [[ ! $1 =~ ^- ]] && min_repeat_length=$1;;
    --max-region-length ) shift; [[ ! $1 =~ ^- ]] && max_region_length=$1;;
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

zgrep --invert-match ^"#" ${repeats}.bed.gz |
  python utilities/get_regular_chromosomes.py |
  python filter-repeats/classify_tandem_repeats_by_length.py ${min_repeat_length} |
  sort --version-sort -k1,1 -k2,2 | # bedtools merge requires sorted input
  bedtools merge -i stdin -c 4 -o collapse | 
  python filter-repeats/classify_merged_tandem_repeats.py |
  awk --assign OFS='\t' '$NF == "1" { print $1, $2, $3 }' | 
  awk --assign max_region_length=${max_region_length} '$3-$2 < max_region_length' 
