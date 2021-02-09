#!/usr/bin/env bash
#SBATCH --time=1:00:00
#SBATCH --cpus-per-task=16
#SBATCH --mem=40g # sacct -o reqmem,maxrss,averss,elapsed -j JOBID
#SBATCH --account=ucgd-rw
#SBATCH --partition=ucgd-rw

experiment=$1
root=$2

reference="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/reference/GRCh38_full_analysis_set_plus_decoy_hla"

population="CHS"
sample="HG00514"

number_threads="16"

bash ${root}/evaluate-calls/evaluate.mantaFlavors.sh \
    --output ${experiment} \
    --threads ${number_threads} \
    --reference ${reference} \
    --population ${population} \
    --sample ${sample} \
    --root ${root} \
  2> ${experiment}/evaluate-calls-mantaFlavors.log
