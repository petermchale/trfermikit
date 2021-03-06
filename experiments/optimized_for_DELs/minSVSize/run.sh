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

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"

genome_build="hg38" # or "hg19"

for minSVSize in 1000 500 250 100 70 60 50; do 
  output="$PWD/data/minSVSize=${minSVSize}"
  mkdir --parents ${output}

  cp ${root}/config.core.json ${output}/config.json
  ${root}/utilities/update_config.sh ${root} ${output} filterCalls minSVSize ${minSVSize}

  ln -s ${root}/experiments/repeats.${genome_build}.tab.gz ${output} 

  sbatch \
    --job-name="minSVSize=${minSVSize}" \
    --output="${output}/slurm.%j.log" \
    run_trfermikit_and_evaluate_calls.sh ${root} ${output} 
done
