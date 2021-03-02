for regions in all genes exons_UTRs; do
  if [[ ${regions} == "all" ]]; then 
    regions="INS"
  fi
  output="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit/experiments/${regions}/singleBaseMatchReward_singleBaseMismatchPenalty_gapOpenPenalties_gapExtensionPenalties/data/gapExtensionPenalties=1,0_gapOpenPenalties=16,41_singleBaseMatchReward=10_singleBaseMismatchPenalty=12"
  bash compute_sv_lengths_of_truth_on_regions.sh --output ${output}
done 