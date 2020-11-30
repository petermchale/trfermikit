#!/usr/bin/env bash

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

population="CHS"
sample="HG00514"

# build 38 of human reference genome:
# *.bed.gz, *.bed.gz.tbi :
repeats="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/repeats/simple-repeats_no-annotations" 
# *.genome, *.mmi, *.fa, *.fa.fai :
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"
# *.cram, *.cram.crai : 
alignments="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/illumina_crams/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/${population}/${sample}/high_cov_alignment/${sample}.alt_bwamem_GRCh38DH.20150715.${population}.high_coverage"
# *.bed.gz, *.bed.gz.tbi :
functional_regions="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/genes/Homo_sapiens.GRCh38.99" 

min_repeat_length="100"

number_threads="16"
svtype="DEL"

output="data/${sample}.${svtype}"

# to facilitate a small-scale test of tool correctness and tool usage: 
repeats_small="${output}/repeats.small"
set +o pipefail
zgrep --invert-match ^"#" ${repeats}.bed.gz | head -1000 | bash utilities/sort_compress_index_regions.sh ${repeats_small}
repeats=${repeats_small}
set -o pipefail

# no need to export PATH since it is already in the environment: https://unix.stackexchange.com/a/26059/406037
PATH="$PWD:$PATH"

# the arguments --min-repeat-length and --functional-regions are optional 
# only svtype==DEL is currently supported
trfermikit \
  --output ${output} \
  --repeats ${repeats} \
  --reference ${reference} \
  --threads ${number_threads} \
  --svtype ${svtype} \
  --alignments ${alignments} \
  --min-repeat-length ${min_repeat_length} \
  --functional-regions ${functional_regions} 

bash evaluate-calls/evaluate.sh \
    --output ${output} \
    --threads ${number_threads} \
    --reference ${reference} \
    --population ${population} \
    --sample ${sample} \
    --svtype ${svtype} \
  2> ${output}/evaluate.log 

