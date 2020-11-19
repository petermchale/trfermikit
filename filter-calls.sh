#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --tr_fermikit_output ) shift; [[ ! $1 =~ ^- ]] && tr_fermikit_output=$1;;
    --svtype ) shift; [[ ! $1 =~ ^- ]] && svtype=$1;;
    --prefix ) shift; [[ ! $1 =~ ^- ]] && prefix=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    --number_threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
    *) bash utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
  esac 
  shift
done

set -o errexit
set -o pipefail
set -o nounset
set -o noclobber

set -o xtrace
# Must use single quote to prevent variable expansion.
# For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# instead of its value when PS4 is used later in the script
# https://stackoverflow.com/a/6697845/6674256
# ${FOO:+val}    val if $FOO is set
# ${FOO[0]}   element #0 of the FOO array
# https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

calls="${tr_fermikit_output}/${prefix}.raw"
unitigs="${tr_fermikit_output}/${prefix}.unsrt"
regions="${tr_fermikit_output}/regions"
slop=500

cat ${calls}.vcf | bash utilities/sort_compress_index_calls.sh ${calls}

# /usr/bin/time won't work with bash functions, only bash scripts
/usr/bin/time --verbose bash utilities/sort_compress_index_alignments.sh ${unitigs}

calls_decomposed_normalized_svtype="${calls}.decomposed.normalized.${svtype}"
bash utilities/decompose_normalize_findSVs.sh \
    --svtype ${svtype} \
    --calls ${calls} \
    --reference ${reference} \
    --number_threads ${number_threads} \
  | bash utilities/sort_compress_index_calls.sh ${calls_decomposed_normalized_svtype}

calls_unitigSupport="${calls_decomposed_normalized_svtype}.unitigSupport"
python filters/filterByUnitigSupport_annotate.py \
    --alignments ${unitigs} \
    --regions ${regions} \
    --calls ${calls_decomposed_normalized_svtype} \
  | bash utilities/sort_compress_index_calls.sh "${calls_unitigSupport}"

calls_thinned="${calls_unitigSupport}.thinned"
bash filters/sparsify_clusters.sh --calls ${calls_unitigSupport} --slop ${slop} \
  | bash utilities/sort_compress_index_calls.sh "${calls_thinned}"
