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

genome_build="hg38" # or "hg19"

# *.fa :
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"
# *.cram, *.cram.crai : 
alignments="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/illumina_crams/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/${population}/${sample}/high_cov_alignment/${sample}.alt_bwamem_GRCh38DH.20150715.${population}.high_coverage"
# *.bed.gz :
functional_regions="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/genes/Homo_sapiens.GRCh38.99" 

min_repeat_length="100"
min_repeat_period="6"

number_threads="16"
svtype="DEL"

trfermikit_path="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
output="${trfermikit_path}/data/${sample}.${svtype}"
PATH="${trfermikit_path}:$PATH"

# the arguments --min-repeat-length, --min-repeat-period, and --functional-regions are optional 
# only svtype==DEL is currently supported
trfermikit \
  --genome-build ${genome_build} \
  --output ${output} \
  --reference ${reference} \
  --threads ${number_threads} \
  --svtype ${svtype} \
  --alignments ${alignments} \
  --min-repeat-length ${min_repeat_length} \
  --min-repeat-period ${min_repeat_period} \
  --functional-regions ${functional_regions} 

bash ${trfermikit_path}/evaluate-calls/evaluate.sh \
    --output ${output} \
    --threads ${number_threads} \
    --reference ${reference} \
    --population ${population} \
    --sample ${sample} \
    --svtype ${svtype} \
    --root ${trfermikit_path} \
  2> ${output}/evaluate-calls.log 
