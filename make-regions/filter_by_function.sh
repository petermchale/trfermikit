#!/usr/bin/env bash

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

functional_regions="$1" 

bash utilities/info.sh "$(bedtools --version)"

# https://github.com/arq5x/bedtools2/issues/834
bedtools intersect -a stdin -b <(zgrep --invert-match ^"#" ${functional_regions}.bed.gz) -wa -u -f 1