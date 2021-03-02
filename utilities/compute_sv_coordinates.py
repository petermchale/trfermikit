from cyvcf2 import VCF
import argparse 
from sv import coordinates

def compute_sv_coordinates():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--calls', type=str, help='')
  args = parser.parse_args()

  variants = VCF('/dev/stdin') if args.calls == 'stdin' else VCF(args.calls)
  
  for variant in variants: 
    chromosome = variant.CHROM    
    start, end = coordinates(variant)
    print('{}\t{}\t{}\n'.format(chromosome, start, end), end='')

  variants.close()

if __name__ == '__main__': 
  compute_sv_coordinates()

