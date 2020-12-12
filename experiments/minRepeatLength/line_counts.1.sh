# set -o xtrace 
# set -o errexit 
# set -o nounset

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
PATH="${root}:$PATH" 

report () {
  local intervals=$1
  local label=$2
  for minRepeatLength in 0 50 100; do
    local path_to_intervals="data/minRepeatLength=${minRepeatLength}/$intervals"
    echo "" 
    echo "$label: minRepeatLength=$minRepeatLength: lines=$(zgrep -v ^# $path_to_intervals | wc -l )"
  done
}
 
report regions.bed.gz "regions"
report fermikit.raw.decomposed.normalized.DEL.vcf.gz "trfermikit calls: all"
report fermikit.raw.decomposed.normalized.DEL.unitigSupport.vcf.gz "trfermikit calls: unitigSupport"
report fermikit.raw.decomposed.normalized.DEL.unitigSupport.thinned.vcf.gz "trfermikit calls: unitigSupport.thinned" 

population="CHS" 
sample="HG00514" 
svtype="DEL"
manta_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/manta/standard_run/results/${population}/${sample}/results/variants/diploidSV"
manta_calls_decomposed_normalized_svtype="${manta_calls}.decomposed.normalized.${svtype}"
echo ""
echo "manta calls: lines=$(zgrep -v ^# ${manta_calls_decomposed_normalized_svtype}.vcf.gz | wc -l )"

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

for minRepeatLength in 0 50 100; do
  regions="data/minRepeatLength=${minRepeatLength}/regions"
  trfermikit_calls="data/minRepeatLength=${minRepeatLength}/fermikit.raw.decomposed.normalized.DEL.unitigSupport.thinned" 
  echo ""
  echo "manta calls in regions: minRepeatLength=$minRepeatLength: lines=$(bedtools intersect -u -f 1 -a $manta_calls_decomposed_normalized_svtype.vcf.gz -b $regions.bed.gz | wc -l)"
  echo "manta calls in pacbio-covered regions: minRepeatLength=$minRepeatLength: lines=$(bedtools intersect -u -f 1 -a $manta_calls_decomposed_normalized_svtype.vcf.gz -b <(pacbio_covered_regions ${regions}) | wc -l)"
  echo "trfermikit calls in regions: minRepeatLength=$minRepeatLength: lines=$(bedtools intersect -u -f 1 -a $trfermikit_calls.vcf.gz -b $regions.bed.gz | wc -l)"
  echo "trfermikit calls in pacbio-covered regions: minRepeatLength=$minRepeatLength: lines=$(bedtools intersect -u -f 1 -a $trfermikit_calls.vcf.gz -b <(pacbio_covered_regions ${regions}) | wc -l)"
  echo "pacbio-covered regions: minRepeatLength=$minRepeatLength: lines=$(pacbio_covered_regions $regions | wc -l)"  
  echo "pacbio calls in pacbio-covered regions: minRepeatLength=$minRepeatLength: lines=$(bedtools intersect -u -f 1 -a $pacbio_calls_decomposed_normalized_svtype.vcf.gz -b <(pacbio_covered_regions ${regions}) | wc -l)"
done

