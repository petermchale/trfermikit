#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --svtype ) shift; [[ ! $1 =~ ^- ]] && svtype=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    --threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
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

calls="${output}/fermikit.raw"
unitigs="${output}/fermikit.srt"
regions="${output}/regions"

parameters=${output}/filter-calls
# Chaisson defines an SV to be an event >50bp in size
# That is, only events >50bp are recorded in the pacbio callset
# Thus, discovered events <50bp may be flagged as FPs by truvari
${root}/bin/jq \
  --null-input \
  '{ 
    "intra cluster distance threshold": "500",
    "minimum SV size": "50",
    "block length threshold": "25",
    "mapping quality threshold": "0"
  }' \
  > ${parameters}.json

calls_decomposed_normalized_svtype="${calls}.decomposed.normalized.${svtype}"
bash ${root}/filter-calls/decompose_normalize_findSVs.sh \
    --svtype ${svtype} \
    --calls ${calls} \
    --reference ${reference} \
    --threads ${number_threads} \
    --parameters ${parameters} \
    --root ${root} \
  | bash ${root}/utilities/sort_compress_index_calls.sh \
    --calls ${calls_decomposed_normalized_svtype} \
    --root ${root}

calls_unitigSupport="${calls_decomposed_normalized_svtype}.unitigSupport"
python ${root}/filter-calls/filterByUnitigSupport_annotate.py \
    --alignments ${unitigs} \
    --regions ${regions} \
    --calls ${calls_decomposed_normalized_svtype} \
    --parameters ${parameters} \
  | bash ${root}/utilities/sort_compress_index_calls.sh \
    --calls "${calls_unitigSupport}" \
    --root ${root}

calls_thinned="${calls_unitigSupport}.thinned"
bash ${root}/filter-calls/sparsify_clusters.sh \
    --calls ${calls_unitigSupport} \
    --parameters ${parameters} \
    --root ${root} \
  | bash ${root}/utilities/sort_compress_index_calls.sh \
    --calls "${calls_thinned}" \
    --root ${root}
