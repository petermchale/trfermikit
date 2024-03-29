#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --calls ) shift; [[ ! $1 =~ ^- ]] && calls=$1;;
    --output ) shift; [[ ! $1 =~ ^- ]] && output=$1;;
    --root ) shift; [[ ! $1 =~ ^- ]] && root=$1;;
    *) echo -e "${RED}$0: $1 is an invalid flag${NO_COLOR}" >&2; exit 1;;
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

read_config () { 
  local key1_=$1 
  local key2_=$2
  ${root}/utilities/read_config.sh ${root} ${output} ${key1_} ${key2_}
}

max_intra_cluster_distance=$(read_config filterCalls maxIntraClusterDistance)

sparsify_clusters () {
  local calls_=$1
  local cluster_column=11
  local confidence_column=12
  ${root}/bin/bedtools cluster -i ${calls_}.vcf.gz -d ${max_intra_cluster_distance} \
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
