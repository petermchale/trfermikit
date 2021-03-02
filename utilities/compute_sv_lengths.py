from cyvcf2 import VCF
import argparse 
from sv import get_sv_length

def compute_sv_lengths():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--calls', type=str, help='')
  args = parser.parse_args()

  variants = VCF('/dev/stdin') if args.calls == 'stdin' else VCF(args.calls)
  
  for variant in variants: 
    print(abs(get_sv_length(variant)))
 
  variants.close()

if __name__ == '__main__': 
  compute_sv_lengths()

