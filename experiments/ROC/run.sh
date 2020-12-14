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

job_count=0 

for minCoverage in 0 5 10; do # trfermikit is "0"
  d1="minCoverage=${minCoverage}"
  for gapOpenPenalties in 3,15 6,26 10,35; do # asm10 is "16,41"; trfermikit is "6,26"
    d2="gapOpenPenalties=${gapOpenPenalties}"
    for singleBaseMatchReward in 10; do # asm10 is "1"; trfermikit is "10"
      # d3="singleBaseMatchReward=${singleBaseMatchReward}"
      for singleBaseMismatchPenalty in 12; do # asm10 is "9"; trfermikit is "12"
        # d4="singleBaseMismatchPenalty=${singleBaseMismatchPenalty}"
        for minUnitigMappingQuality in 0 10; do # trfermikit is "0" 
          d5="minUnitigMappingQuality=${minUnitigMappingQuality}"
          for minUnitigBlockLength in 15 25 40; do # trfermikit is "25" 
            d6="minUnitigBlockLength=${minUnitigBlockLength}"
            experiment="${d1};${d2};${d5};${d6}"
            output="${root}/experiments/ROC/data/${experiment}"
            mkdir --parents ${output}
            cp ${root}/config.core.json ${output}/config.json
            ${root}/utilities/update_config.sh ${root} ${output} makeRegions minCoverage ${minCoverage}
            ${root}/utilities/update_config.sh ${root} ${output} makeCalls gapOpenPenalties ${gapOpenPenalties}
            ${root}/utilities/update_config.sh ${root} ${output} makeCalls singleBaseMatchReward ${singleBaseMatchReward}
            ${root}/utilities/update_config.sh ${root} ${output} makeCalls singleBaseMismatchPenalty ${singleBaseMismatchPenalty}
            ${root}/utilities/update_config.sh ${root} ${output} filterCalls minUnitigMappingQuality ${minUnitigMappingQuality}
            ${root}/utilities/update_config.sh ${root} ${output} filterCalls minUnitigBlockLength ${minUnitigBlockLength}
            
            sbatch \
              --job-name="${experiment}" \
              --output="${output}/slurm.%j.log" \
              ${root}/experiments/ROC/run_trfermikit_and_evaluate_calls.sh ${root} ${output}
            ((job_count++))          
          done
        done
      done
    done
  done
done
 
bash ${root}/utilities/info.sh "number of jobs submitted: $job_count"
