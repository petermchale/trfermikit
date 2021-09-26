#!/usr/bin/env bash
#SBATCH --time=6:00:00
#SBATCH --nodes=1
#SBATCH --mem=40g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID
# a slurm task is a Linux process:
# https://support.ceci-hpc.be/doc/_contents/QuickStart/SubmittingJobs/SlurmTutorial.html#going-parallel
#SBATCH --ntasks=16
# slurm does not allocate resources for more than 16 CPUs per job,
# so request one CPU per Linux process:
#SBATCH --cpus-per-task=1 
#SBATCH --account=quinlan-rw
#SBATCH --partition=quinlan-shared-rw

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

root=$1
output=$2 
population_sample="${3}"

PATH="${root}:$PATH"

reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"

IFS=, read population sample <<< "${population_sample}"

alignments="/scratch/ucgd/lustre-work/quinlan/data-shared/HGSVC2/illumina-alignments/${sample}.final"

number_threads="16"

minRepeatLength="100" 

trfermikit \
  --output ${output} \
  --reference ${reference} \
  --threads ${number_threads} \
  --alignments ${alignments} \
  --min-repeat-length ${minRepeatLength} \

bash ${root}/evaluate-calls/evaluate.stratify_by_region_length.HGSVC2.sh \
    --output ${output} \
    --population ${population} \
    --sample ${sample} \
    --root ${root} \
  2> ${output}/evaluate.stratify_by_region_length.log 

