#!/usr/bin/env bash

export CYAN='\033[0;36m'
export RED='\033[0;31m'
export NO_COLOR='\033[0m'

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --population ) shift; [[ ! $1 =~ ^- ]] && population=$1;;
    --sample ) shift; [[ ! $1 =~ ^- ]] && sample=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    *) echo -e "${RED}$0: $1 is an invalid flag${NO_COLOR}" >&2; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber

# set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

evaluate_calls () {
  local svtype_=$1
  local region_size_range_="${2}"
  bash ${root}/evaluate-calls/evaluate.stratify_by_region_length.HGSVC2.core.sh \
      --output ${output} \
      --population ${population} \
      --sample ${sample} \
      --root ${root} \
      --svtype ${svtype_} \
      --region-size-range ${region_size_range_} \
    2> ${output}/evaluate.stratify_by_region_length.${svtype_}.${region_size_range_}.log
}

read_config () { 
  local key1_=$1 
  local key2_=$2
  ${root}/utilities/read_config.sh ${root} ${output} ${key1_} ${key2_}
}

slop=$(read_config makeRegions slop)
bash ${root}/utilities/info.sh "shortest region is $((2 * $slop))"
bash ${root}/utilities/info.sh "longest region is $(read_config makeRegions maxRegionLength)"

for svtype in "DEL" "INS"; do 
  for region_size_range in "600,625" "625,650" "650,700" "700,800" "800,1000" "1000,2000" "2000,100000"; do 
    evaluate_calls ${svtype} ${region_size_range}
  done
done 

