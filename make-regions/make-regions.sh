#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
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

read_config () { 
  local key1_=$1 
  local key2_=$2
  ${root}/utilities/read_config.sh ${root} ${output} ${key1_} ${key2_}
}

functional_regions=$(read_config makeRegions functionalRegions)
genome_build=$(read_config makeRegions genomeBuild)

reference=$(read_config general reference)
number_threads=$(read_config general numberThreads)
alignments=$(read_config general alignments) 

slop=$(read_config makeRegions slop)
min_coverage=$(read_config makeRegions minCoverage)
max_coverage=$(read_config makeRegions maxCoverage)

if [[ "${functional_regions}" == "none" ]]; then
  overlapped_functional_regions=false
else
  overlapped_functional_regions=true
fi
${root}/utilities/update_config.sh ${root} ${output} makeRegions overlappedFunctionalRegions ${overlapped_functional_regions}

repeats="${output}/repeats.${genome_build}"
if [[ ! -f ${repeats}.tab.gz ]]; then 
  bash ${root}/make-regions/download_simple_repeats.sh \
    --root ${root} \
    --output ${output} \
    --repeats ${repeats}
else 
  bash ${root}/utilities/info.sh "skipping the downloading of repeats"
fi 

filter_repeats_by_length_and_function () {
  bash ${root}/make-regions/filter_by_length.sh \
    --repeats ${repeats} \
    --output ${output} \
    --root ${root} \
  | 
  # https://unix.stackexchange.com/a/38311/406037 : 
  if [[ "${functional_regions}" == "none" ]]; then 
    cat
  else 
    bash ${root}/make-regions/filter_by_function.sh \
      --root ${root} \
      --output ${output}
  fi 
}

scratch=$(mktemp --tmpdir=${output} --directory)
clean_up () {
  rm --recursive --force "${scratch}"
}
trap clean_up EXIT

mosdepth_prefix="${scratch}/mosdepth.coverage"
${root}/bin/mosdepth \
  --no-per-base \
  --fast-mode \
  --threads ${number_threads} \
  --by <(filter_repeats_by_length_and_function) \
  --fasta ${reference}.fa \
  ${mosdepth_prefix} \
  ${alignments}.cram

${root}/bin/bedtools slop -i ${mosdepth_prefix}.regions.bed.gz -g ${reference}.genome -b ${slop} \
  | awk --assign min_coverage=${min_coverage} '$4 > min_coverage' \
  | awk --assign max_coverage=${max_coverage} '$4 < max_coverage' \
  | awk --assign OFS='\t' '{ print $1, $2, $3 }' \
  | bash ${root}/utilities/sort_compress_index_regions.sh \
    --regions "${output}/regions" \
    --root ${root}


