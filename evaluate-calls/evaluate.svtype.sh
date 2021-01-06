#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    --population ) shift; [[ ! $1 =~ ^- ]] && population=$1;;
    --sample ) shift; [[ ! $1 =~ ^- ]] && sample=$1;;
    --svtype ) shift; [[ ! $1 =~ ^- ]] && svtype=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    *) bash ${root}/utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
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

#######################################################

# https://stackoverflow.com/a/43476575/6674256
export PYTHONPATH="${root}/utilities"

#######################################################

truvari_trfermikit="truvari-${svtype}-pacbio-trfermikit"
truvari_manta="truvari-${svtype}-pacbio-manta"

pacbio_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/calls/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/${sample}.BIP-unified.filtered"
tr_fermikit_calls="${output}/fermikit.raw.decomposed.normalized.${svtype}"
manta_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/manta/standard_run/results/${population}/${sample}/results/variants/diploidSV"

pacbio_calls_decomposed_normalized_svtype="${output}/pacbioCalls.decomposed.normalized.${svtype}"
manta_calls_decomposed_normalized_svtype="${output}/mantaCalls.decomposed.normalized.${svtype}"

parameters="${output}/config"

#######################################################

bash ${root}/filter-calls/decompose_normalize_findSVs.sh \
    --svtype ${svtype} \
    --calls ${pacbio_calls} \
    --reference ${reference} \
    --threads ${number_threads} \
    --parameters ${parameters} \
    --root ${root} \
  | bash ${root}/utilities/sort_compress_index_calls.sh \
    --calls ${pacbio_calls_decomposed_normalized_svtype} \
    --root ${root}

bash ${root}/filter-calls/decompose_normalize_findSVs.sh \
    --svtype ${svtype} \
    --calls ${manta_calls} \
    --reference ${reference} \
    --threads ${number_threads} \
    --parameters ${parameters} \
    --root ${root} \
  | bash ${root}/utilities/sort_compress_index_calls.sh \
    --calls ${manta_calls_decomposed_normalized_svtype} \
    --root ${root}

pacbio_covered_regions () { 
  local pacbio_covered_regions_on_h0_="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/pacbio_local_assemblies/${sample}.h0.covered.sorted"
  local pacbio_covered_regions_on_h1_="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/pacbio_local_assemblies/${sample}.h1.covered.sorted"

  local regions_="${output}/regions"

  ${root}/bin/bedtools intersect -a ${regions_}.bed.gz -b ${pacbio_covered_regions_on_h0_}.bed -wa -u -f 1 |
    ${root}/bin/bedtools intersect -a stdin -b ${pacbio_covered_regions_on_h1_}.bed -wa -u -f 1 |
    sort --version-sort -k1,1 -k2,2
} 

run_truvari () { 
  local comp_calls_="$1"
  local truvari_output_="$2"
  
  rm --recursive --force ${truvari_output_}

  # comp file argument must be bgzipped: 
  # https://github.com/spiralgenetics/truvari/blob/e1d9a4ea441b8102fb3e352ac6fcbf65fe9405e2/truvari/truvari#L705
  # this version of truvari resolved the following issue: 
  # https://github.com/spiralgenetics/truvari/issues/49
  /usr/bin/time --verbose truvari bench \
    --prog \
    --debug \
    --sizemin 30 \
    --pctsize 0.25 \
    --reference ${reference}.fa \
    --base ${pacbio_calls_decomposed_normalized_svtype}.vcf.gz \
    --comp ${comp_calls_}.vcf.gz \
    --output ${truvari_output_} \
    --includebed <(pacbio_covered_regions) \
    > ${truvari_output_}.log 2>&1
} 

run_truvari \
  "${tr_fermikit_calls}" \
  "${output}/${truvari_trfermikit}"
 
run_truvari \
  "${tr_fermikit_calls}.unitigSupport" \
  "${output}/${truvari_trfermikit}.unitigSupport"

run_truvari \
  "${tr_fermikit_calls}.unitigSupport.thinned" \
  "${output}/${truvari_trfermikit}.unitigSupport.thinned"

run_truvari \
  "${manta_calls_decomposed_normalized_svtype}" \
  "${output}/${truvari_manta}"
