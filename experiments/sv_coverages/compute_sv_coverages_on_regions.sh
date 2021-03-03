#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --calls-name ) shift; [[ ! $1 =~ ^- ]] && calls_name=$1;;
    --svtype ) shift; [[ ! $1 =~ ^- ]] && svtype=$1;;
    *) echo -e "${RED}$0: $1 is an invalid flag${NO_COLOR}" >&2; exit 1;;
  esac 
  shift
done

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

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"

# https://stackoverflow.com/a/43476575/6674256
export PYTHONPATH="${root}/utilities"

population="CHS"
sample="HG00514"
alignments="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/illumina_crams/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/${population}/${sample}/high_cov_alignment/${sample}.alt_bwamem_GRCh38DH.20150715.${population}.high_coverage"

if [[ ${calls_name} == "pacbio" ]]; then 
  calls="${output}/pacbioCalls.decomposed.normalized.${svtype}"
elif [[ ${calls_name} == "trfermikit" ]]; then 
  calls="${output}/fermikit.raw.decomposed.normalized.${svtype}.unitigSupport.thinned"
elif [[ ${calls_name} == "trfermikit_TP" ]]; then 
  calls="${output}/truvari-${svtype}-pacbio-trfermikit.unitigSupport.thinned/tp-call"
elif [[ ${calls_name} == "trfermikit_FP" ]]; then 
  calls="${output}/truvari-${svtype}-pacbio-trfermikit.unitigSupport.thinned/fp"
elif [[ ${calls_name} == "trfermikit_FN" ]]; then 
  calls="${output}/truvari-${svtype}-pacbio-trfermikit.unitigSupport.thinned/fn"
elif [[ ${calls_name} == "manta" ]]; then 
  calls="${output}/mantaCalls.decomposed.normalized.${svtype}"
else
  echo -e "${RED}${calls_name} is invalid${NO_COLOR}" >&2
fi

pacbio_covered_regions () { 
  local pacbio_covered_regions_on_h0_="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/pacbio_local_assemblies/${sample}.h0.covered.sorted"
  local pacbio_covered_regions_on_h1_="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/pacbio_local_assemblies/${sample}.h1.covered.sorted"

  local regions_="${output}/regions"

  ${root}/bin/bedtools intersect -a ${regions_}.bed.gz -b ${pacbio_covered_regions_on_h0_}.bed -wa -u -f 1 |
    ${root}/bin/bedtools intersect -a stdin -b ${pacbio_covered_regions_on_h1_}.bed -wa -u -f 1 |
    sort --version-sort -k1,1 -k2,2
} 

calls_file () {
  if [[ -f ${calls}.vcf.gz ]]; then 
    echo ${calls}.vcf.gz
  elif [[ -f ${calls}.vcf ]]; then
    echo ${calls}.vcf
  else 
    echo -e "${RED}Neither ${calls}.vcf.gz nor ${calls}.vcf exists!${NO_COLOR}" >&2
  fi 
}

sv_regions () {
  ${root}/bin/bedtools intersect \
      -a $(calls_file) \
      -b <(pacbio_covered_regions) \
      -wa -u -f 1 -header \
    | python ${root}/utilities/compute_sv_coordinates.py --calls stdin \
    | ${root}/bin/bedtools slop -i stdin -g ${reference}.genome -b 200
}

scratch=$(mktemp --tmpdir=${PWD} --directory)
clean_up () {
  rm --recursive --force "${scratch}"
}
trap clean_up EXIT

mosdepth_prefix="${scratch}/mosdepth.coverage"
${root}/bin/mosdepth \
  --no-per-base \
  --fast-mode \
  --threads 16 \
  --by <(sv_regions) \
  --fasta ${reference}.fa \
  ${mosdepth_prefix} \
  ${alignments}.cram

zcat ${mosdepth_prefix}.regions.bed.gz 




