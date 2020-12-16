root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
svtype="DEL"
output="minRepeatLength/data/minRepeatLength=0" 
parameters="${output}/config.json"
tr_fermikit_calls="${output}/fermikit.raw"

# https://stackoverflow.com/a/43476575/6674256
export PYTHONPATH="${root}/utilities"

zcat ${tr_fermikit_calls}.vcf.gz | 
  python ${root}/filter-calls/find_SVs.py \
      --calls stdin \
      --svtype ${svtype} \
      --parameters ${parameters}