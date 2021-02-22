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

export CYAN='\033[0;36m'
export RED='\033[0;31m'
export NO_COLOR='\033[0m'

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"

job_count=0 

for gapExtensionPenalties in 1,0 2,1 3,2; do 
  d1="gapExtensionPenalties=${gapExtensionPenalties}"
  for gapOpenPenalties in 6,26 16,41 25,50; do 
    d2="gapOpenPenalties=${gapOpenPenalties}"
    for singleBaseMatchReward in 1 10 20; do 
      d3="singleBaseMatchReward=${singleBaseMatchReward}"
      for singleBaseMismatchPenalty in 5 9 12; do 
        d4="singleBaseMismatchPenalty=${singleBaseMismatchPenalty}"
        experiment="${d1}_${d2}_${d3}_${d4}"
        output="$PWD/data/${experiment}"
        mkdir --parents ${output}

        cp ${root}/config.json ${output}/config.json
        ${root}/utilities/update_config.sh ${root} ${output} makeCalls gapExtensionPenalties ${gapExtensionPenalties}
        ${root}/utilities/update_config.sh ${root} ${output} makeCalls gapOpenPenalties ${gapOpenPenalties}
        ${root}/utilities/update_config.sh ${root} ${output} makeCalls singleBaseMatchReward ${singleBaseMatchReward}
        ${root}/utilities/update_config.sh ${root} ${output} makeCalls singleBaseMismatchPenalty ${singleBaseMismatchPenalty}

        ln -s ${root}/experiments/repeats.hg38.tab.gz ${output} 

        sbatch \
          --job-name="${experiment}" \
          --output="${output}/slurm.%j.log" \
          $PWD/run_trfermikit_and_evaluate_calls.sh ${root} ${output}

        ((job_count++))          
      done
    done
  done
done

bash ${root}/utilities/info.sh "number of jobs submitted: $job_count"