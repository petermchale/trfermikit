set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

population="CHS" 
sample="HG00514" 

svtype="DEL" 

output="minCoverage_gapOpenPenalties_minUnitigMappingQuality_minUnitigBlockLength/data/minCoverage=0_gapOpenPenalties=5,20_minUnitigMappingQuality=10_minUnitigBlockLength=25"

pacbio_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/calls/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/${sample}.BIP-unified.filtered"
tr_fermikit_calls="${output}/fermikit.raw"
manta_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/manta/standard_run/results/${population}/${sample}/results/variants/diploidSV"

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"
number_threads="16" 

# https://stackoverflow.com/a/43476575/6674256
export PYTHONPATH="${root}/utilities"

compute_min_SV_size () {
  local calls_=$1 

  ${root}/bin/bcftools norm \
      --check-ref w \
      --fasta-ref ${reference}.fa \
      --multiallelics \
      -any \
      --threads ${number_threads} \
      ${calls_}.vcf.gz \
      2> /dev/null 
    # | python compute_min_SV_size.py \
    #   --svtype ${svtype} \
    #   --calls "stdin"
}

echo "min-size pacbio call of type ${svtype}:"
compute_min_SV_size ${pacbio_calls}
echo "min-size trfermikit call of type ${svtype}:"
compute_min_SV_size ${tr_fermikit_calls}
echo "min-size manta call of type ${svtype}:"
compute_min_SV_size ${manta_calls}

