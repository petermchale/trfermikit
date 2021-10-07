import sys 
from cyvcf2 import VCF
import json 

def count_TP_FP_FN(directory_combined_caller, type_combined_caller, directory_individual_caller):
  vcf = VCF(f'{directory_combined_caller}/{type_combined_caller}.sorted.vcf.gz')

  TP = 0
  FP = 0
  for variant in vcf: 
    if variant.INFO.get('TruScore'): 
      TP += 1
    else: 
      FP += 1

  vcf.close()

  with open(f'{directory_individual_caller}/summary.txt', 'r') as fh_individual: 
    counts = json.load(fh_individual)
    event_count = counts['TP-base'] + counts['FN']

  with open(f'{directory_combined_caller}/counts.json', 'w') as fh_combined: 
    json.dump({
      'TP-base': TP, # evaluate.ipynb assumes the existence of this key 
      'FP': FP,
      'FN': event_count - TP , # evaluate.ipynb assumes the existence of this key 
    }, fh_combined, indent=2)

if __name__ == '__main__':
  count_TP_FP_FN(
    directory_combined_caller=sys.argv[1],
    type_combined_caller=sys.argv[2],
    directory_individual_caller=sys.argv[3]
  )