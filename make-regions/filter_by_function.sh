#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    *) echo -e "${RED}$0: $1 is an invalid flag${NO_COLOR}" >&2; exit 1;;
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

functional_regions=$(read_config makeRegions functionalRegions)
minRepeatFunctionalOverlap=$(read_config makeRegions minRepeatFunctionalOverlap)

bash ${root}/utilities/info.sh "$(${root}/bin/bedtools --version)"

# https://github.com/arq5x/bedtools2/issues/834
${root}/bin/bedtools intersect -a stdin -b <(zgrep --invert-match ^"#" ${functional_regions}.bed.gz) -wa -u -f ${minRepeatFunctionalOverlap}