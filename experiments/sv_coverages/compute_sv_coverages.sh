#!/usr/bin/env bash
#SBATCH --time=0:30:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID
#SBATCH --account=ucgd-rw
#SBATCH --partition=ucgd-rw

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --alignments-name ) shift; [[ ! $1 =~ ^- ]] && alignments_name=$1;;
    *) echo -e "${RED}$0: $1 is an invalid flag${NO_COLOR}" >&2; exit 1;;
  esac 
  shift
done

export CYAN='\033[0;36m'
export RED='\033[0;31m'
export NO_COLOR='\033[0m'

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

svtype="DEL"

for calls_name in trfermikit_TP trfermikit_FP trfermikit_FN; do
  for regions in regions_intersecting_exons_and_UTRs regions_intersecting_genes all_regions; do
    if [[ ${regions} == "all_regions" ]]; then 
      directory="INS"
    elif [[ ${regions} == "regions_intersecting_genes" ]]; then
      directory="genes"
    elif [[ ${regions} == "regions_intersecting_exons_and_UTRs" ]]; then
      directory="exons_UTRs"
    else
      echo -e "${RED}${regions} is invalid${NO_COLOR}" >&2
    fi
    output="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit/experiments/${directory}/singleBaseMatchReward_singleBaseMismatchPenalty_gapOpenPenalties_gapExtensionPenalties/data/gapExtensionPenalties=1,0_gapOpenPenalties=16,41_singleBaseMatchReward=10_singleBaseMismatchPenalty=12"
    bash compute_sv_coverages_on_regions.sh \
        --output ${output} \
        --svtype ${svtype} \
        --calls-name ${calls_name} \
        --alignments-name ${alignments_name} \
      > sv_coverages.${regions}.${svtype}.${calls_name}.${alignments_name}.csv
  done
done 