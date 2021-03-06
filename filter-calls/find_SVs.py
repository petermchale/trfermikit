from cyvcf2 import VCF
import sys 
import argparse 
from sv import get_sv_length, hom_ref, get_svtype

def find_SVs():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--calls', type=str, help='')
  parser.add_argument('--svtype', type=str, help='')
  parser.add_argument('--parameters', type=str, help='')
  args = parser.parse_args()

  import json 
  parameters = json.load(open('{}.json'.format(args.parameters)))

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

    minSVSize = int(parameters['filterCalls']['minSVSize'])
    if get_svtype(variant) == svtype and abs(get_sv_length(variant)) >= minSVSize: 
      print(variant, end='')
 
  variants.close()

if __name__ == '__main__': 
  find_SVs()

