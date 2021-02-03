population="CHS"
sample="HG00514"
alignments="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/illumina_crams/ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/hgsv_sv_discovery/data/${population}/${sample}/high_cov_alignment/${sample}.alt_bwamem_GRCh38DH.20150715.${population}.high_coverage"
region="chr1:10000-10100"

fetch_reads_in_region () { 
  samtools view -h ${alignments}.cram ${region}
} 

unmapped_read_in_region=$(fetch_reads_in_region | samtools view -f 4 | cut -f1 | head -1)

echo "an unmapped read associated with given region that is present in fastq file (and therefore not indicated as unmapped):"
fetch_reads_in_region | samtools fastq | grep ${unmapped_read_in_region} 
