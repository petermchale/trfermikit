#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --svtype ) shift; [[ ! $1 =~ ^- ]] && svtype=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    --threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
    *) bash utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
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

cluster_distance="500"

# Chaisson defines an SV to be an event >50bp in size
# That is, only events >50bp are recorded in the pacbio callset
# Thus, discovered events <50bp may be flagged as FPs by truvari
sv_length_threshold="50"

jq \
  --null-input \
  --arg cluster_distance ${cluster_distance} \
  --arg sv_length_threshold ${sv_length_threshold} \
  '{ 
    "intra cluster distance threshold": $cluster_distance,
    "minimum SV size": $sv_length_threshold
  }' \
  > ${output}/filter-calls.json

gunzip --force --stdout ${calls}.vcf.gz | bash utilities/sort_compress_index_calls.sh ${calls}

calls_decomposed_normalized_svtype="${calls}.decomposed.normalized.${svtype}"
bash filter-calls/decompose_normalize_findSVs.sh \
    --svtype ${svtype} \
    --calls ${calls} \
    --reference ${reference} \
    --number_threads ${number_threads} \
    --sv-length-threshold ${sv_length_threshold} \
  | bash utilities/sort_compress_index_calls.sh ${calls_decomposed_normalized_svtype}

calls_unitigSupport="${calls_decomposed_normalized_svtype}.unitigSupport"
python filter-calls/filterByUnitigSupport_annotate.py \
    --alignments ${unitigs} \
    --regions ${regions} \
    --calls ${calls_decomposed_normalized_svtype} \
  | bash utilities/sort_compress_index_calls.sh "${calls_unitigSupport}"

calls_thinned="${calls_unitigSupport}.thinned"
bash filter-calls/sparsify_clusters.sh --calls ${calls_unitigSupport} --cluster-distance ${cluster_distance} \
  | bash utilities/sort_compress_index_calls.sh "${calls_thinned}"
