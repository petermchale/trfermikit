#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --calls ) shift; [[ ! $1 =~ ^- ]] && calls=$1;;
    --parameters ) shift; [[ ! $1 =~ ^- ]] && parameters=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    --intra-cluster-distance-threshold ) shift; [[ ! $1 =~ ^- ]] && intra_cluster_distance_threshold=$1;;
    *) bash ${root}/utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

sparsify_clusters () {
  local calls_=$1
  local cluster_column=11
  local confidence_column=12
  ${root}/bin/bedtools cluster -i ${calls_}.vcf.gz -d ${intra_cluster_distance_threshold} \
    | python ${root}/filter-calls/append_INFO_value_to_vcf_record.py "Confidence" \
    | sort -k${cluster_column},${cluster_column}n -k${confidence_column},${confidence_column}nr \
    | ${root}/bin/bedtools groupby -grp ${cluster_column} -opCols ${confidence_column} -ops max -full \
    | cut -f1-$((cluster_column-1))
}

vcf_headers () {
  local calls_=$1
  zgrep ^"#" ${calls_}.vcf.gz
}

(vcf_headers ${calls}; sparsify_clusters ${calls})
