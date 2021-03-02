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

# https://stackoverflow.com/a/43476575/6674256
export PYTHONPATH="${root}/utilities"

sample="HG00514"

if [[ ${calls_name} == "pacbio" ]]; then 
  calls="${output}/pacbioCalls.decomposed.normalized.${svtype}"
elif [[ ${calls_name} == "trfermikit" ]]; then 
  calls="${output}/fermikit.raw.decomposed.normalized.${svtype}.unitigSupport.thinned"
elif [[ ${calls_name} == "trfermikit_TP" ]]; then 
  calls="${output}/truvari-${svtype}-pacbio-trfermikit.unitigSupport.thinned/tp-call"
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

${root}/bin/bedtools intersect \
    -a $(calls_file) \
    -b <(pacbio_covered_regions) \
    -wa -u -f 1 -header \
  | python ${root}/utilities/compute_sv_lengths.py --calls stdin 

