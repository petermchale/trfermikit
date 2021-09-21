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
# output="${root}/experiments/INS/singleBaseMatchReward_singleBaseMismatchPenalty_gapOpenPenalties_gapExtensionPenalties/data/gapExtensionPenalties=1,0_gapOpenPenalties=16,41_singleBaseMatchReward=1_singleBaseMismatchPenalty=12"
population_sample="${3}"

PATH="${root}:$PATH"

reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"

IFS=, read population sample <<< "${population_sample}"

alignments="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/illumina_crams/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/${population}/${sample}/high_cov_alignment/${sample}.alt_bwamem_GRCh38DH.20150715.${population}.high_coverage"

number_threads="16"

minRepeatLength="100" 

trfermikit \
  --output ${output} \
  --reference ${reference} \
  --threads ${number_threads} \
  --alignments ${alignments} \
  --min-repeat-length ${minRepeatLength} \

bash ${root}/evaluate-calls/evaluate.stratify_by_region_length.sh \
    --output ${output} \
    --population ${population} \
    --sample ${sample} \
    --root ${root} \
  2> ${output}/evaluate.stratify_by_region_length.log 
