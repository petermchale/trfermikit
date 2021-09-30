consortium="HGSVC2"
population="GWD"
sample="HG02818"

output="${consortium}/data/${population},${sample}"

for caller in "manta" "trfermikit.unitigSupport.thinned"; do
  for svtype in "DEL" "INS"; do 
    for region_size_range in "600,625" "625,650" "650,700" "700,800" "800,1000" "1000,2000" "2000,100000"; do 
      truvari="${output}/truvari-${svtype}-${region_size_range}-pacbio-${caller}"

      # TODO: check headers of combined callset 
      bcftools concat ${truvari}/fp.vcf ${truvari}/tp-call.vcf > ${truvari}/combined-callset-${caller}.vcf 
    done
  done 
done


