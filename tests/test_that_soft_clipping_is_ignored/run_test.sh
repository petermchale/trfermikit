#set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber

# set -o xtrace
# # Must use single quote to prevent variable expansion.
# # For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# # instead of its value when PS4 is used later in the script
# # https://stackoverflow.com/a/6697845/6674256
# # ${FOO:+val}    val if $FOO is set
# # ${FOO[0]}   element #0 of the FOO array
# # https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
# PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"
region_with_soft_clips_in_unitigs="chr9:125,985,047-125,985,344"
experiment="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit/experiments/INS/singleBaseMatchReward_singleBaseMismatchPenalty_gapOpenPenalties_gapExtensionPenalties/data/gapExtensionPenalties=1,0_gapOpenPenalties=6,26_singleBaseMatchReward=1_singleBaseMismatchPenalty=12"

${root}/make-calls/fermi.kit/htsbox pileup -c -f ${reference}.fa -r ${region_with_soft_clips_in_unitigs} \
  ${experiment}/fermikit.srt.bam 
