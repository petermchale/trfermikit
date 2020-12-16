population="CHS" 
sample="HG00514" 

svtype="DEL" 

output="minRepeatLength/data/minRepeatLength=0" 

pacbio_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/calls/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/${sample}.BIP-unified.filtered"
tr_fermikit_calls="${output}/fermikit.raw"
manta_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/manta/standard_run/results/${population}/${sample}/results/variants/diploidSV"

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"
number_calls="16" 

compute_min_SV_size () {
  local calls_=$1 

  ${root}/bin/bcftools norm \
      --check-ref w \
      --fasta-ref ${reference}.fa \
      --multiallelics \
      -any \
      --threads ${number_threads} \
      ${calls_}.vcf.gz \
    | python compute_min_SV_size.py ${svtype}
}

echo "min size of pacbio calls of type ${svtype} is: $(compute_min_SV_size ${pacbio_calls})" 
echo "min size of trfermikit calls of type ${svtype} is: $(compute_min_SV_size ${tr_fermikit_calls})" 
echo "min size of manta calls of type ${svtype} is: $(compute_min_SV_size ${manta_calls})" 

