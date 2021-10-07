set -o errexit 
# set -o xtrace 
set -o nounset 

root="/scratch/ucgd/lustre-work/quinlan/u6018199/chaisson_2019/analysis/locally_assemble_short_reads/trfermikit"

overlap_fraction="0.9"
svtype="DEL"

create-callset () { 
  local dir_=$1 

  cat ${dir_}/fp.vcf \
    | ${root}/utilities/sort_compress_index_calls.sh \
      --calls ${dir_}/fp.sorted \
      --root ${root} \
    2> /dev/null 
  cat ${dir_}/tp-call.vcf \
    | ${root}/utilities/sort_compress_index_calls.sh \
      --calls ${dir_}/tp-call.sorted \
      --root ${root} \
    2> /dev/null 
  ${root}/bin/bcftools concat --allow-overlaps \
      ${dir_}/fp.sorted.vcf.gz \
      ${dir_}/tp-call.sorted.vcf.gz \
    2> /dev/null \
    | ${root}/utilities/sort_compress_index_calls.sh \
      --calls ${dir_}/callset.sorted \
      --root ${root} \
    2> /dev/null 
}

manta-intersect-trfermikit () {
  ${root}/bin/bedtools intersect -header -u -f ${overlap_fraction} -r \
      -a ${VNTR_category_manta}/callset.sorted.vcf.gz \
      -b ${VNTR_category_trfermikit}/callset.sorted.vcf.gz \
    | ${root}/utilities/sort_compress_index_calls.sh \
      --calls ${VNTR_category_combined}/manta-intersect-trfermikit.sorted \
      --root ${root} \
    2> /dev/null 
}

manta-less-trfermikit () { 
  $root/bin/bedtools subtract -header -A -f ${overlap_fraction} -r \
      -a ${VNTR_category_manta}/callset.sorted.vcf.gz \
      -b ${VNTR_category_trfermikit}/callset.sorted.vcf.gz \
    | ${root}/utilities/sort_compress_index_calls.sh \
      --calls ${VNTR_category_combined}/manta-less-trfermikit.sorted \
      --root ${root} \
    2> /dev/null 
}

trfermikit-less-manta () { 
  $root/bin/bedtools subtract -header -A -f ${overlap_fraction} -r \
      -a ${VNTR_category_trfermikit}/callset.sorted.vcf.gz \
      -b ${VNTR_category_manta}/callset.sorted.vcf.gz \
    | python relabel-sample.py ${sample} \
    | ${root}/utilities/sort_compress_index_calls.sh \
      --calls ${VNTR_category_combined}/trfermikit-less-manta.sorted \
      --root ${root} \
    2> /dev/null 
}

manta-union-trfermikit () { 
  ${root}/bin/bcftools concat --allow-overlaps \
      ${VNTR_category_combined}/manta-intersect-trfermikit.sorted.vcf.gz \
      ${VNTR_category_combined}/manta-less-trfermikit.sorted.vcf.gz \
      ${VNTR_category_combined}/trfermikit-less-manta.sorted.vcf.gz \
    2> /dev/null \
    | ${root}/utilities/sort_compress_index_calls.sh \
      --calls ${VNTR_category_union}/manta-union-trfermikit.sorted \
      --root ${root} \
    2> /dev/null 

  # TODO: remove 
  rm --force ${VNTR_category_combined}/combined.sorted.vcf.gz*
  rm --force ${VNTR_category_combined}/counts.json
  
}

# use a glob pattern to locate the sample directories:  
for path in */data/*; do
  IFS=, read prefix sample <<< ${path}
  echo "sample: ${sample}"  
  for region_size_range in "600,625" "625,650" "650,700" "700,800" "800,1000" "1000,2000" "2000,100000"; do 
    echo -e "\tregion_size_range: ${region_size_range}"

    VNTR_category_manta="${path}/truvari-${svtype}-${region_size_range}-pacbio-manta"
    VNTR_category_trfermikit="${path}/truvari-${svtype}-${region_size_range}-pacbio-trfermikit.unitigSupport.thinned"
    VNTR_category_combined="${path}/truvari-${svtype}-${region_size_range}-pacbio-combined"
    VNTR_category_union="${path}/truvari-${svtype}-${region_size_range}-pacbio-union"
    VNTR_category_intersection="${path}/truvari-${svtype}-${region_size_range}-pacbio-intersection"

    create-callset ${VNTR_category_manta}
    create-callset ${VNTR_category_trfermikit} 

    mkdir --parents ${VNTR_category_combined}
    manta-intersect-trfermikit
    manta-less-trfermikit
    trfermikit-less-manta

    mkdir --parents ${VNTR_category_union}
    manta-union-trfermikit
    python count-TP-FP-FN.py ${VNTR_category_union} "manta-union-trfermikit" ${VNTR_category_manta}

    mkdir --parents ${VNTR_category_intersection}
    cp ${VNTR_category_combined}/manta-intersect-trfermikit.sorted.vcf.gz* ${VNTR_category_intersection}
    python count-TP-FP-FN.py ${VNTR_category_intersection} "manta-intersect-trfermikit" ${VNTR_category_manta}
  done
done 


