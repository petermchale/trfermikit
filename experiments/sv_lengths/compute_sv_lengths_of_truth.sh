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

for regions in all genes exons_UTRs; do
  if [[ ${regions} == "all" ]]; then 
    directory="INS"
  else
    directory=${regions}
  fi
  output="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit/experiments/${directory}/singleBaseMatchReward_singleBaseMismatchPenalty_gapOpenPenalties_gapExtensionPenalties/data/gapExtensionPenalties=1,0_gapOpenPenalties=16,41_singleBaseMatchReward=10_singleBaseMismatchPenalty=12"
  bash compute_sv_lengths_of_truth_on_regions.sh --output ${output} \
    > sv_lengths.${regions}.txt
done 