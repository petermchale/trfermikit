from cyvcf2 import VCF
import pysam
import sys
import numpy as np 
import argparse
import gzip
from sv import coordinates

def parse(locus): 
  chromosome, start_end = locus.split(':')
  start, end = map(lambda s: int(s.replace(',', '')), start_end.split('-'))
  return chromosome, start, end

def length(block): 
  start, end = block 
  return end - start

def max_block(blocks): 
  block_sizes = [length(block) for block in blocks]  
  largest_block_index = np.argmax(block_sizes)
  return block_sizes[largest_block_index], largest_block_index

def retainCall_reportConfidence(unitigs, variant, region, parameters): 
  min_mapping_quality = int(parameters['filterCalls']['minUnitigMappingQuality'])
  min_unitig_block_length = int(parameters['filterCalls']['minUnitigBlockLength'])

  call_start, call_end = coordinates(variant) 
 
  for unitig in unitigs.fetch(*parse(region)):
    if unitig.mapping_quality <= min_mapping_quality:
      continue
    # https://pysam.readthedocs.io/en/latest/api.html#pysam.AlignedSegment.get_blocks
    blocks = unitig.get_blocks()
    for index in range(1,len(blocks)):
      previous_block = blocks[index-1]
      current_block = blocks[index]
      if previous_block[1] == call_start and current_block[0] == call_end:
        blocks_upstream_of_call = blocks[:index]
        blocks_downstream_of_call = blocks[index:]
        max_block_size_upstream_of_call, max_block_index_upstream_of_call = max_block(blocks_upstream_of_call) 
        max_block_size_downstream_of_call, max_block_index_downstream_of_call = max_block(blocks_downstream_of_call) 
#        block_immediately_upstream_of_call = blocks_upstream_of_call[-1]
#        block_immediately_downstream_of_call = blocks_downstream_of_call[0]
#        condition_1 = length(block_immediately_upstream_of_call) > min_unitig_block_length
#        condition_2 = length(block_immediately_downstream_of_call) > min_unitig_block_length
        condition_1 = max_block_size_upstream_of_call > min_unitig_block_length
        condition_2 = max_block_size_downstream_of_call > min_unitig_block_length
        # condition_3 = max_block_index_upstream_of_call == 0
        # condition_4 = max_block_index_downstream_of_call == len(blocks_downstream_of_call) - 1
        call_confidence = (max_block_size_upstream_of_call + max_block_size_downstream_of_call)/len(blocks)
        # Adding condition_3 and condition_4 to condition_1 and condition_2 drastically reduces FP count while minimally impacting TP count.
        # This suggests that a RNN might give better overall performance than the current mash-up of heuristic filters.
        # if condition_1 and condition_2 and condition_3 and condition_4: 
        # if condition_3 and condition_4: 
        if condition_1 and condition_2: 
          return True, call_confidence 

  return False, None

def annotate(variant, confidence): 
  INFO_index = 7
  line = str(variant)
  fields = line.split('\t')
  INFO = fields[INFO_index]
  key_value = '{}={}'.format('Confidence',confidence)
  if INFO == '.':
    INFO = key_value
  else:
    INFO += ';' + key_value
  fields[INFO_index] = INFO
  return '\t'.join(fields)

def filter_annotate_calls(): 
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--alignments', type=str, help='')
  parser.add_argument('--regions', type=str, help='')
  parser.add_argument('--calls', type=str, help='')
  parser.add_argument('--parameters', type=str, help='')
  args = parser.parse_args()

  import json 
  parameters = json.load(open('{}.json'.format(args.parameters)))

  vcf = VCF(args.calls+'.vcf.gz')
  vcf.add_info_to_header({
    'ID': 'Confidence', 
    'Description': 'Measure of confidence in call based upon unitig structure',
    'Type': 'String', 
    'Number': '1'
  })
  with pysam.AlignmentFile(args.alignments+'.bam', 'rb') as unitigs, gzip.open(args.regions+'.bed.gz', 'rt') as regions:
    print(vcf.raw_header, end='')
    for region in regions:
      chromosome, start, end = region.strip().split('\t') 
      region = '{}:{}-{}'.format(chromosome, start, end)
      for variant in vcf(region): 
        retain_call, call_confidence = retainCall_reportConfidence(unitigs, variant, region, parameters)
        if retain_call:
          print(annotate(variant, call_confidence), end='') 
  vcf.close()

if __name__ == "__main__":
  filter_annotate_calls()
