set -o errexit
set -o pipefail
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

export RED='\033[0;31m'
export CYAN='\033[0;36m'
export NO_COLOR='\033[0m'

run_test () {
  local test_=$1
  local test_directory_=$PWD/test_${test_}
  bash ${test_directory_}/test_trfermikit.sh ${test_directory_} \
    > ${test_directory_}/expected_vs_observed.txt \
    2> ${test_directory_}/test_trfermikit.log
}

run_test calls 
run_test functional_regions 

bash test_that_soft_clipping_is_ignored/run_test.sh > test_that_soft_clipping_is_ignored/expected_vs_observed.txt 2>&1

cat test_calls/expected_vs_observed.txt \
    test_functional_regions/expected_vs_observed.txt \
    test_that_soft_clipping_is_ignored/expected_vs_observed.txt \
  | less -R 
