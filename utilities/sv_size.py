from color_text import error, info
import sys 

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

