# set -o xtrace 
# set -o errexit 
# set -o nounset

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
PATH="${root}:$PATH" 

population="CHS" 
sample="HG00514" 
svtype="DEL"

manta_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/manta/standard_run/results/${population}/${sample}/results/variants/diploidSV"
manta_calls_decomposed_normalized_svtype="${manta_calls}.decomposed.normalized.${svtype}"

pacbio_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/calls/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/${sample}.BIP-unified.filtered"
pacbio_calls_decomposed_normalized_svtype="${pacbio_calls}.decomposed.normalized.${svtype}"

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
  regions="data/minRepeatLength=${minRepeatLength}/regions"
  trfermikit_calls="data/minRepeatLength=${minRepeatLength}/fermikit.raw.decomposed.normalized.DEL.unitigSupport.thinned" 
  echo ""
  echo "minRepeatLength=$minRepeatLength"
  echo "--------------------------------"
  echo "# manta calls: $(number_calls $manta_calls_decomposed_normalized_svtype $regions)" 
  echo "# trfermikit calls: $(number_calls $trfermikit_calls $regions)" 
  echo "# pacbio calls: $(number_calls $pacbio_calls_decomposed_normalized_svtype $regions)" 
  echo "# regions: $(pacbio_covered_regions ${regions} | wc -l)"
done

