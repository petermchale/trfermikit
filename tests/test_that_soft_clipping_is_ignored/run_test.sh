#set -o errexit
set -o pipefail
set -o nounset
# set -o noclobber

# set -o xtrace
# # Must use single quote to prevent variable expansion.
# # For example, if double quotes were used, ${LINENO} would take on the value of the current line,
# # instead of its value when PS4 is used later in the script
# # https://stackoverflow.com/a/6697845/6674256
# # ${FOO:+val}    val if $FOO is set
# # ${FOO[0]}   element #0 of the FOO array
# # https://www.gnu.org/software/bash/manual/html_node/Bash-Variables.html
# PS4='+ (${BASH_SOURCE[0]##*/} @ ${LINENO}): ${FUNCNAME[0]:+${FUNCNAME[0]}(): }'

export CYAN='\033[0;36m'
export RED='\033[0;31m'
export NO_COLOR='\033[0m'

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"
region_with_soft_clips_in_unitigs="chr9:125,985,047-125,985,344"

experiment="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit/experiments/INS/singleBaseMatchReward_singleBaseMismatchPenalty_gapOpenPenalties_gapExtensionPenalties/data/gapExtensionPenalties=1,0_gapOpenPenalties=6,26_singleBaseMatchReward=1_singleBaseMismatchPenalty=12"
unitigs="${experiment}/fermikit.srt"

population="CHS"
sample="HG00514"
pacbio_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/calls/ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/${sample}.BIP-unified.filtered"
manta_calls="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/manta/standard_run/results/${population}/${sample}/results/variants/diploidSV"

echo "***************************"
echo -e "${CYAN}test_that_soft_clipping_is_ignored: ${NO_COLOR}\n"

bash ${root}/utilities/info.sh "htsbox-pileup calls in region where unitigs are soft-clipped:"
${root}/make-calls/fermi.kit/htsbox pileup -c -f ${reference}.fa -r ${region_with_soft_clips_in_unitigs} ${unitigs}.bam \
  | grep -v ^# 

bash ${root}/utilities/info.sh "pacbio INSs in region:"
zgrep 125985180 ${pacbio_calls}.vcf.gz
zgrep 125985189 ${pacbio_calls}.vcf.gz

bash ${root}/utilities/info.sh "manta INSs in region:"
zgrep 125985168 ${manta_calls}.vcf.gz
