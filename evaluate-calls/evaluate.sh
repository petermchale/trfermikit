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

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

evaluate_calls () {
  local svtype=$1
  bash ${root}/evaluate-calls/evaluate.svtype.sh \
      --output ${output} \
      --population ${population} \
      --sample ${sample} \
      --svtype ${svtype} \
      --root ${root} \
    2> ${output}/evaluate-calls-${svtype}.log
}

evaluate_calls "DEL"
evaluate_calls "INS"
