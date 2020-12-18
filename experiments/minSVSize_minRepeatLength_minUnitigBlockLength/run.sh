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

genome_build="hg38" # or "hg19"

job_count=0 

for minSVSize in 500 100 50; do 
  d1="minSVSize=${minSVSize}"
  for minRepeatLength in 100 50 0; do 
    d2="minRepeatLength=${minRepeatLength}"
    for minUnitigBlockLength in 15 25 40; do 
      d3="minUnitigBlockLength=${minUnitigBlockLength}"

      experiment="${d1}_${d2}_${d3}"
      output="${root}/experiments/minSVSize_minRepeatLength_minUnitigBlockLength/data/${experiment}"
      mkdir --parents ${output}

      cp ${root}/config.core.json ${output}/config.json
      ${root}/utilities/update_config.sh ${root} ${output} filterCalls minSVSize ${minSVSize}
      ${root}/utilities/update_config.sh ${root} ${output} filterCalls minUnitigBlockLength ${minUnitigBlockLength}

      ln -s ${root}/experiments/repeats.${genome_build}.tab.gz ${output} 

      sbatch \
        --job-name="${experiment}" \
        --output="${output}/slurm.%j.log" \
        ${root}/experiments/minSVSize_minRepeatLength_minUnitigBlockLength/run_trfermikit_and_evaluate_calls.sh ${root} ${output} ${genome_build} ${minRepeatLength}
      ((job_count++))          
    done
  done
done
 
bash ${root}/utilities/info.sh "number of jobs submitted: $job_count"
