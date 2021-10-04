import sys 
from cyvcf2 import VCF

def print_vcf(sample):
  vcf = VCF('/dev/stdin')

  for header_line in vcf.raw_header.split('\n'):
    if header_line.startswith('#CHROM'): continue
    if len(header_line) == 0 : continue
    print(header_line)    
  print('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}'.format(sample))

  for variant in vcf: print(variant, end='')  
  
  vcf.close()
  
if __name__ == '__main__':
  print_vcf(sample=sys.argv[1])