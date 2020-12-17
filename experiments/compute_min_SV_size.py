from cyvcf2 import VCF
import sys 
from sv import get_sv_length, hom_ref, get_svtype
import argparse 

def compute_min_SV_size():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--calls', type=str, help='')
  parser.add_argument('--svtype', type=str, help='')
  args = parser.parse_args()

  variants = VCF('/dev/stdin') if args.calls == 'stdin' else VCF(args.calls)
  
  svtype = args.svtype
  if svtype not in ['DEL', 'INS']:
    print('svtype', svtype, 'not permitted!', file=sys.stderr) 
    sys.exit(1) 

  min_size_variant = None
  min_size = 10000

  for variant in variants: 
    # Decomposition may cause vcf records with genotype "0/0" to appear. 
    # These should be removed 
    if hom_ref(variant): 
      continue

    if get_svtype(variant) == svtype and abs(get_sv_length(variant)) < min_size:
      min_size_variant = variant 
      min_size = get_sv_length(variant)
      
  variants.close()

  print('min-size variant: {}'.format(str(min_size_variant)), end='')
  print('min size: {}'.format(min_size))
  print()

if __name__ == '__main__': 
  compute_min_SV_size()