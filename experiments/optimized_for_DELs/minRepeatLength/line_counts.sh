# set -o xtrace 
# set -o errexit 
# set -o nounset

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
PATH="${root}:$PATH" 

svtype=$1 
# echo "svtype: ${svtype}" 
# echo "--------------------------------"

population="CHS" 
sample="HG00514" 

pacbio_covered_regions () { 
  local pacbio_covered_regions_on_h0_="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/pacbio_local_assemblies/${sample}.h0.covered.sorted"
  local pacbio_covered_regions_on_h1_="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/pacbio_local_assemblies/${sample}.h1.covered.sorted"

  local regions_=$1

  bedtools intersect -a ${regions_}.bed.gz -b ${pacbio_covered_regions_on_h0_}.bed -wa -u -f 1 |
    bedtools intersect -a stdin -b ${pacbio_covered_regions_on_h1_}.bed -wa -u -f 1 |
    sort --version-sort -k1,1 -k2,2
} 

number_calls () { 
  local calls_=$1
  local regions_=$2
  bedtools intersect -u -f 1 -a ${calls_}.vcf.gz -b <(pacbio_covered_regions ${regions_}) \
    | wc -l
} 

for minRepeatLength in 0 50 100; do
  experiment="data/minRepeatLength=${minRepeatLength}"
  regions="${experiment}/regions"
  trfermikit_calls="${experiment}/fermikit.raw.decomposed.normalized.${svtype}.unitigSupport.thinned" 
  manta_calls="${experiment}/mantaCalls.decomposed.normalized.${svtype}"
  pacbio_calls="${experiment}/pacbioCalls.decomposed.normalized.${svtype}"
  echo ""
  echo "minRepeatLength=$minRepeatLength"
  echo "--------------------------------"
  echo "# manta calls: $(number_calls $manta_calls $regions)" 
  echo "# trfermikit calls: $(number_calls $trfermikit_calls $regions)" 
  echo "# pacbio calls: $(number_calls $pacbio_calls $regions)" 
  echo "# regions: $(pacbio_covered_regions ${regions} | wc -l)"
done

