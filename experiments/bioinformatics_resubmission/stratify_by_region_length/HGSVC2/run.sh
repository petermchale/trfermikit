set -o errexit 
set -o nounset 

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"

jq --raw-output '.[]' "/scratch/ucgd/lustre-work/quinlan/data-shared/HGSVC2/calls/evaluation-samples.json" \
  | while read population_sample; do
    output="$PWD/data/${population_sample}"
    mkdir --parents ${output}

    if [[ ! -f ${output}/repeats.hg38.tab.gz ]]; then
      ln -s ${root}/experiments/repeats.hg38.tab.gz ${output} 
    fi

    sbatch \
      --job-name="HGSVC2.stratify_by_region_length.${population_sample}" \
      --output="${output}/slurm.log" \
      $PWD/trfermikit.evaluate.sh ${root} ${output} ${population_sample}
  done
