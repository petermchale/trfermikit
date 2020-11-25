#!/usr/bin/env bash

# https://devhints.io/bash#miscellaneous
# put option-fetching before "set -o nounset" so that we can detect flags without arguments
while [[ "$1" =~ ^- ]]; do 
  case $1 in
    --svtype ) shift; [[ ! $1 =~ ^- ]] && svtype=$1;;
    --calls ) shift; [[ ! $1 =~ ^- ]] && calls=$1;;
    --reference ) shift; [[ ! $1 =~ ^- ]] && reference=$1;;
    --number_threads ) shift; [[ ! $1 =~ ^- ]] && number_threads=$1;;
    --sv-length-threshold ) shift; [[ ! $1 =~ ^- ]] && sv_length_threshold=$1;;
    *) bash utilities/error.sh "$0: $1 is an invalid flag"; exit 1;;
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

# There will be N ALTs if there are N overlapping high-base-quality unitigs/reads, 
# each with different gaps, say.
# For a diploid sample, however, at most two of these ALTs may appear in the genotype. 
# So split ("decompose") multiallelic sites into separate biallelic vcf records, and trim and left-align variants:
# https://genome.sph.umich.edu/wiki/Variant_Normalization
# "normalization will only be applied if the --fasta-ref option is supplied" to "bcftools norm":
# http://samtools.github.io/bcftools/bcftools.html
# also see: http://www.cureffi.org/2014/04/24/converting-genetic-variants-to-their-minimal-representation/
bcftools norm --check-ref w --fasta-ref ${reference}.fa --multiallelics -any --threads ${number_threads} ${calls}.vcf.gz \
  | python filter-calls/find_SVs.py \
    --calls stdin \
    --svtype ${svtype} \
    --sv-length-threshold ${sv_length_threshold}



