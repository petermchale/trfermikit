from cyvcf2 import VCF
import sys 
from color_text import error, info
import color_traceback
import argparse 

def hom_ref(variant): 
  # http://brentp.github.io/cyvcf2/#cyvcf2
  if len(variant.genotypes) > 1: 
    error('multi-sample vcf record:') 
    info(str(variant)) 
    sys.exit(1) 
  genotype, = variant.genotypes
  allele_haplotype_0, allele_haplotype_1, _ = genotype
  return True if allele_haplotype_0 == 0 and allele_haplotype_1 == 0 else False

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

def get_svtype(variant): 
  if variant.INFO.get('SVTYPE'): 
    # pacbio and manta vcf:
    return variant.INFO.get('SVTYPE')
  else: 
    # tr-fermikit vcf: 
    return 'DEL' if get_sv_length(variant) < 0 else 'INS'

def find_SVs():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--calls', type=str, help='')
  parser.add_argument('--svtype', type=str, help='')
  # https://docs.python.org/3/library/argparse.html#dest : 
  parser.add_argument('--sv-length-threshold', dest='sv_length_threshold', type=int, help='')
  args = parser.parse_args()

  variants = VCF('/dev/stdin') if args.calls == 'stdin' else VCF(args.calls)
  
  svtype = args.svtype
  if svtype not in ['DEL', 'INS']:
    print('svtype', svtype, 'not permitted!', file=sys.stderr) 
    sys.exit(1) 
 
  print(variants.raw_header, end="")  
  for variant in variants: 
    # Decomposition may cause vcf records with genotype "0/0" to appear. 
    # These should be removed because truvari flags these as FPs, artificially inflating the FP rate
    if hom_ref(variant): 
      continue

    if get_svtype(variant) == svtype and abs(get_sv_length(variant)) >= args.sv_length_threshold: 
      print(variant, end='')
 
  variants.close()

if __name__ == '__main__': 
  find_SVs()

