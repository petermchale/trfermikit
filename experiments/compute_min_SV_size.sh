population="CHS" 
sample="HG00514" 

svtype="DEL" 

output="minRepeatLength/data/minRepeatLength=0" 

pacbio_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/calls/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/${sample}.BIP-unified.filtered"
tr_fermikit_calls="${output}/fermikit.raw"
manta_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/manta/standard_run/results/${population}/${sample}/results/variants/diploidSV"

echo "min size of pacbio calls of type ${svtype} is: $(python compute_min_SV_size.py ${pacbio_calls}.vcf.gz)" 
echo "min size of trfermikit calls of type ${svtype} is: $(python compute_min_SV_size.py ${tr_fermikit_calls}.vcf.gz)" 
echo "min size of manta calls of type ${svtype} is: $(python compute_min_SV_size.py ${manta_calls}.vcf.gz)" 

