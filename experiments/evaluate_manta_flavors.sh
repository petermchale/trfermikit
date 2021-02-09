root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"

experiments () {
  ls -d \
    INS/singleBaseMatchReward_singleBaseMismatchPenalty_gapOpenPenalties_gapExtensionPenalties/data/* \
    optimized_for_DELs/minCoverage_gapOpenPenalties_minUnitigMappingQuality_minUnitigBlockLength/data/*
} 

job_count=0

for experiment in $(experiments); do 
#  sbatch \
#    --job-name="${experiment}" \
#    --output="${output}/slurm.mantaFlavors.%j.log" \
#    evaluate_manta_flavors_core.sh ${experiment} ${root} 

  evaluate_manta_flavors_core.sh ${experiment} ${root}

  ((job_count++))
done

bash ${root}/utilities/info.sh "number of jobs submitted: ${job_count}"
