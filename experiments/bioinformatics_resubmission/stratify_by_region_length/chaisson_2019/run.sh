set -o errexit 
set -o nounset 

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"

for population_sample in \
  "CHS,HG00514" \
  "PUR,HG00733" \
  "YRI,NA19240" ; do
    output="$PWD/data/${population_sample}"
    mkdir --parents ${output}

    if [[ ! -f ${output}/repeats.hg38.tab.gz ]]; then
      ln -s ${root}/experiments/repeats.hg38.tab.gz ${output} 
    fi

    sbatch \
      --job-name="stratify_by_region_length.${population_sample}" \
      --output="${output}/slurm.log" \
      $PWD/trfermikit.evaluate.sh ${root} ${output} ${population_sample}
done

