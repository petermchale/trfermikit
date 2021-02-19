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

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"
PATH="${root}:$PATH"

output="results"
reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"
number_threads="16"
alignments="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit/tests/alignments" 
functional_regions="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/genes/Homo_sapiens.GRCh38.99"

rm -rf ${output} 
mkdir --parents ${output}

repeats="repeats.hg38"
if [[ -e ${repeats}.tab.gz ]]; then
  ln -s ${PWD}/${repeats}.tab.gz ${output} 
fi

trfermikit \
  --output ${output} \
  --reference ${reference} \
  --threads ${number_threads} \
  --alignments ${alignments} \
  --functional-regions ${functional_regions}

regions="${output}/regions"

expected_regions='chr1\t853165\t854156\n'
observed_regions=$(zcat ${regions}.bed.gz) 

echo -e "expected regions:\n${expected_regions}" > expected_vs_observed.txt
echo -e "observed regions:\n${observed_regions}" >> expected_vs_observed.txt

cat expected_vs_observed.txt
