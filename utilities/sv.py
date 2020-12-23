from color_text import error, info
import sys 
import color_traceback

def get_sv_length(variant): 
  if variant.INFO.get('SVLEN'): 
    # pacbio and manta vcf:
    return variant.INFO.get('SVLEN') 
  else: 
    # tr-fermikit vcf: 
    if len(variant.ALT) > 1:
      error('There is more than one ALT allele!')
      error('Please decompose the variant:')
      info(str(variant))
      sys.exit(1)
    return len(variant.ALT[0]) - len(variant.REF)

def hom_ref(variant): 
  # http://brentp.github.io/cyvcf2/#cyvcf2
  if len(variant.genotypes) > 1: 
    error('multi-sample vcf record:') 
    info(str(variant)) 
    sys.exit(1) 
  genotype, = variant.genotypes
  allele_haplotype_0, allele_haplotype_1, _ = genotype
  return True if allele_haplotype_0 == 0 and allele_haplotype_1 == 0 else False

def get_svtype(variant): 
  if variant.INFO.get('SVTYPE'): 
    # pacbio and manta vcf:
    return variant.INFO.get('SVTYPE')
  else: 
    # tr-fermikit vcf: 
    if get_sv_length(variant) < 0: return 'DEL' 
    if get_sv_length(variant) > 0: return 'INS'
    if get_sv_length(variant) == 0: return 'SNP'


