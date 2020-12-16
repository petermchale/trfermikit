from cyvcf2 import VCF
import sys 

for variant in VCF(sys.argv[1]): 

sample = 'HG00514'
experiment = 5
number_calls = 0
for variant in VCF('nstd152.GRCh38.variant_call.vcf.gz'):
  if variant.INFO.get('EXPERIMENT') == experiment and variant.INFO.get('SAMPLE') == sample:
    #print(variant)
    number_calls += 1
print('number of calls where experiment = {} in {}: {}'.format(experiment, sample, number_calls))

union = 'PacBio,Illumina'
sv_type = '<DEL>'
number_calls = 0
for variant in VCF('ftp.ncbi.nlm.nih.gov/pub/dbVar/data/Homo_sapiens/by_study/genotype/nstd152/' + sample + '.BIP-unified.vcf.gz'):
  if variant.ALT[0] == sv_type and variant.INFO.get('UNION') == union: 
    #print(variant)
    number_calls += 1
print('number of calls where union = {} and sv_type = {} in {}: {}'.format(union, sv_type, sample, number_calls))

