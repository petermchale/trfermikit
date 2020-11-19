#!/usr/bin/env bash

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
prefix="trfermikit"

tr_fermikit_output="data/${sample}.${svtype}"

# to facilitate a small-scale test of tool correctness and tool usage: 
repeats_small="${tr_fermikit_output}/repeats.small"
set +o pipefail
zgrep --invert-match ^"#" ${repeats}.bed.gz | head -100 | bash utilities/sort_compress_index_regions.sh ${repeats_small}
repeats=${repeats_small}
set -o pipefail

export PATH="$PWD:$PATH"

# the arguments --min-repeat-length and --functional-regions are optional 
# only svtype==DEL is currently supported
trfermikit \
  --trfermikit-output ${tr_fermikit_output} \
  --repeats ${repeats} \
  --reference ${reference} \
  --number-threads ${number_threads} \
  --svtype ${svtype} \
  --prefix ${prefix} \
  --alignments ${alignments} \
  --min-repeat-length ${min_repeat_length} 
# --functional-regions ${functional_regions} 

# evaluate-trfermikit \
#   --trfermikit-output ${tr_fermikit_output} \
#   --number-threads ${number_threads} \
#   --reference ${reference} \
#   --population ${population} \
#   --sample ${sample} \
#   --svtype ${svtype} \
#   --prefix ${prefix}