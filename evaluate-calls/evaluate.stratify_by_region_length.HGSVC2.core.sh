#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --population ) shift; [[ ! $1 =~ ^- ]] && population=$1;;
    --sample ) shift; [[ ! $1 =~ ^- ]] && sample=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    --svtype ) shift; [[ ! $1 =~ ^- ]] && svtype=$1;;
    --region-size-range ) shift; [[ ! $1 =~ ^- ]] && region_size_range=$1;;
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

#######################################################

read_config () { 
  local key1_=$1 
  local key2_=$2
  ${root}/utilities/read_config.sh ${root} ${output} ${key1_} ${key2_}
}

reference=$(read_config general reference)

#######################################################

# https://stackoverflow.com/a/43476575/6674256
export PYTHONPATH="${root}/utilities"

#######################################################

truvari_trfermikit="truvari-${svtype}-${region_size_range}-pacbio-trfermikit"
truvari_manta="truvari-${svtype}-${region_size_range}-pacbio-manta"

pacbio_calls="/scratch/ucgd/lustre-work/quinlan/data-shared/HGSVC2/calls/variants_freeze4_sv_insdel_alt.${population}.${sample}"
tr_fermikit_calls="${output}/fermikit.raw.decomposed.normalized.${svtype}"
manta_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/manta/standard_run/results_bioinformatics_resubmission/${population}/${sample}/results/variants/diploidSV"

pacbio_calls_decomposed_normalized_svtype="${output}/pacbioCalls.decomposed.normalized.${svtype}"
manta_calls_decomposed_normalized_svtype="${output}/mantaCalls.decomposed.normalized.${svtype}"

parameters="${output}/config"

IFS=, read min_region_length max_region_length <<< "${region_size_range}"

#######################################################

bash ${root}/filter-calls/decompose_normalize_findSVs.sh \
    --svtype ${svtype} \
    --calls ${pacbio_calls} \
    --parameters ${parameters} \
    --root ${root} \
    --output ${output} \
  | bash ${root}/utilities/sort_compress_index_calls.sh \
    --calls ${pacbio_calls_decomposed_normalized_svtype} \
    --root ${root}

bash ${root}/filter-calls/decompose_normalize_findSVs.sh \
    --svtype ${svtype} \
    --calls ${manta_calls} \
    --parameters ${parameters} \
    --root ${root} \
    --output ${output} \
  | bash ${root}/utilities/sort_compress_index_calls.sh \
    --calls ${manta_calls_decomposed_normalized_svtype} \
    --root ${root}

pacbio_covered_regions () { 
  local pacbio_covered_regions_on_h1_="/scratch/ucgd/lustre-work/quinlan/data-shared/HGSVC2/hifi-covered-regions/${sample}.hifi.h1.covered.sorted"
  local pacbio_covered_regions_on_h2_="/scratch/ucgd/lustre-work/quinlan/data-shared/HGSVC2/hifi-covered-regions/${sample}.hifi.h2.covered.sorted"

  local regions_="${output}/regions"

  ${root}/bin/bedtools intersect -a ${regions_}.bed.gz -b ${pacbio_covered_regions_on_h1_}.bed -wa -u -f 1 |
    ${root}/bin/bedtools intersect -a stdin -b ${pacbio_covered_regions_on_h2_}.bed -wa -u -f 1 |
    sort --version-sort -k1,1 -k2,2
} 

y_chromosome () {
  echo -e "chrY\t1\t57227415"
}

white_list () {
  # assumes stdin (and stdout) is "chr <TAB> start <TAB> end"
  low_confidence_regions="/scratch/ucgd/lustre-work/quinlan/data-shared/HGSVC2/exclude-regions/LowConfidenceRegions.sorted"
  ${root}/bin/bedtools subtract -a stdin -b ${low_confidence_regions}.bed -A \
    | ${root}/bin/bedtools subtract -a stdin -b <(y_chromosome) -A \
    | sort --version-sort -k1,1 -k2,2
}

# echo "test the white_list function: "
# test_regions () {
#   echo -e "chr1\t121609189\t121609289" 
#   echo -e "chr1\t121609789\t121609889"
#   echo -e "chrY\t1\t100"
# }
# echo "test regions: "
# test_regions 
# echo "test regions that are white listed: "
# test_regions | white_list 
# exit 1 

filter_regions_by_length () { 
  pacbio_covered_regions \
    | white_list \
    | awk --assign max_region_length=${max_region_length} '$3-$2 < max_region_length' \
    | awk --assign min_region_length=${min_region_length} '$3-$2 > min_region_length'
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
    --includebed <(filter_regions_by_length) \
    > ${truvari_output_}.log 2>&1

  filter_regions_by_length \
    | bash ${root}/utilities/sort_compress_index_regions.sh \
      --regions "${truvari_output_}/regions" \
      --root ${root}
} 

run_truvari \
  "${tr_fermikit_calls}.unitigSupport.thinned" \
  "${output}/${truvari_trfermikit}.unitigSupport.thinned"

run_truvari \
  "${manta_calls_decomposed_normalized_svtype}" \
  "${output}/${truvari_manta}"


