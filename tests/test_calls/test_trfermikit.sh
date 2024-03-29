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

test_directory=$1 

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
PATH="${root}:$PATH"

output="${test_directory}/results"
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"
number_threads="16"
alignments="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit/tests/alignments" 

rm -rf ${output} 
mkdir --parents ${output}

trfermikit \
  --output ${output} \
  --reference ${reference} \
  --threads ${number_threads} \
  --alignments ${alignments} 

calls="${output}/fermikit.raw.decomposed.normalized.DEL.unitigSupport.thinned" 

expected_call="chr1    820928  .       CTTGTCCAGCAGGTCCACCCTGTCTACACTACCTGCCTGCAAAGCAGATCCACCCTGTCTACACTACCTGG C       43      .       Confidence=26.25        GT:AD   1/0:14,43" 
observed_call=$(zgrep -v ^# ${calls}.vcf.gz)

echo "***************************" 
echo -e "${CYAN}test_calls: ${NO_COLOR}" 
echo -e "expected call:\n${expected_call}" 
echo -e "observed call:\n${observed_call}" 

