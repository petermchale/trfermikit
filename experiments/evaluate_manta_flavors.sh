# set -o errexit
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

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"

experiments () {
  ls -d \
    INS/singleBaseMatchReward_singleBaseMismatchPenalty_gapOpenPenalties_gapExtensionPenalties/data/* \
    optimized_for_DELs/minCoverage_gapOpenPenalties_minUnitigMappingQuality_minUnitigBlockLength/data/*
} 

job_count=0

for experiment in $(experiments); do 
  sbatch \
    --job-name="${experiment}" \
    --output="${experiment}/slurm.mantaFlavors.%j.log" \
    evaluate_manta_flavors_core.sh ${experiment} ${root} 

  # bash evaluate_manta_flavors_core.sh ${experiment} ${root} 

  ((job_count++))
done

bash ${root}/utilities/info.sh "number of jobs submitted: ${job_count}"
