#!/usr/bin/env bash
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID
#SBATCH --account=ucgd-rw
#SBATCH --partition=ucgd-rw

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

experiment=$1
root=$2

reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"

population="CHS"
sample="HG00514"

number_threads="16"

bash ${root}/evaluate-calls/evaluate.mantaFlavors.sh \
    --output ${experiment} \
    --threads ${number_threads} \
    --reference ${reference} \
    --population ${population} \
    --sample ${sample} \
    --root ${root} \
  2> ${experiment}/evaluate-calls-mantaFlavors.log
