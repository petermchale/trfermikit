#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --repeats ) shift; [[ ! $1 =~ ^- ]] && repeats=$1;;
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

read_config () { 
  local key1_=$1 
  local key2_=$2
  ${root}/utilities/read_config.sh ${root} ${output} ${key1_} ${key2_}
}

min_repeat_length=$(read_config makeRegions minRepeatLength)
max_region_length=$(read_config makeRegions maxRegionLength)
min_repeat_period=$(read_config makeRegions minRepeatPeriod)

exit 1

zgrep --invert-match ^"#" ${repeats}.tab.gz |
  python ${root}/make-regions/classify_tandem_repeats_by_length.py \
    --min-repeat-length ${min_repeat_length} \
    --min-repeat-period ${min_repeat_period} \
    --log-file ${repeats}.classified.tab.gz |
  python ${root}/utilities/get_regular_chromosomes.py |
  sort --version-sort -k1,1 -k2,2 | # bedtools merge requires sorted input
  ${root}/bin/bedtools merge -i stdin -c 4 -o collapse | # assumes that column 4 contains classification
  python ${root}/make-regions/classify_merged_tandem_repeats.py |
  awk --assign OFS='\t' '$NF == "1" { print $1, $2, $3 }' | 
  awk --assign max_region_length=${max_region_length} '$3-$2 < max_region_length' 
