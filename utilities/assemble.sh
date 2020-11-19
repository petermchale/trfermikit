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

EXE_FERMI2=$1
K_UNITIG=$2 
K_MERGE=$3 
N_THREADS=$4 
dependency=$5
target=$6 

# these functions anticipate the error: 
# fermi.pre.gz.log: fermi2: unitig.c:414: fm6_unitig: Assertion `e->mcnt[1] >= n_threads * 2' failed.
# return status of function is return status of last executed command in function

turn_off_parallelization () {
  grep --quiet "Turn off parallelization for this batch as too few strings are left." "${dependency}.log"
}

filtered_fastq_empty () {
  # Returns true when fermi.flt.fq.gz is empty 
  # The file fermi.flt.fq.gz contains “trimmed” reads 
  # obtained from the file fermi.ec.fq.gz by running “bfc -1”, c.f., 
  # https://github.com/lh3/bfc/tree/a73dad248dc56d9d4d22eacbbbc51ac276045168#usage

  grep --quiet "\[M::main_ropebwt2\] symbol counts: ($, A, C, G, T, N) = (0, 0, 0, 0, 0, 0)" "${dependency}.log"
  # NOTE: [ and ] are special characters to grep
}

rm --force filtered_fastq_empty
if filtered_fastq_empty; then 
  touch filtered_fastq_empty
  exit 1 # this causes make to exit 
fi 

if turn_off_parallelization; then 
  number_threads=1
else 
  number_threads=${N_THREADS}
fi

/usr/bin/time --verbose ${EXE_FERMI2} assemble -l ${K_UNITIG} -m ${K_MERGE} -t ${number_threads} ${dependency} 2> "${target}.log" | gzip -1 > ${target} 


