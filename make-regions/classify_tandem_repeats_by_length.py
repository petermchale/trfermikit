import sys
import color_traceback
import argparse
from color_text import error

def classify(fields, args): 
  _, start, end, period = fields
  condition1 = int(end) - int(start) > args.min_repeat_length
  condition2 = int(period) > args.min_repeat_period
  return 1 if condition1 and condition2 else 0

def add_column():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--min-repeat-length', type=int, dest='min_repeat_length', help='')
  parser.add_argument('--min-repeat-period', type=int, dest='min_repeat_period', help='')  
  args = parser.parse_args()

  for line_number, line in enumerate(sys.stdin):
    fields = line.strip().split()
    print(line_number, fields)
    if line_number == 0:
      condition1 = fields[0] != 'chrom' 
      condition2 = fields[1] != 'chromStart'
      condition3 = fields[2] != 'chromEnd'
      condition4 = fields[3] != 'period'      
      if condition1 or condition2 or condition3 or condition4: 
        error('repeats file does not have correct column names')
        sys.exit(1)
    else: 
      new_fields = fields + [classify(fields, args)]
      print('\t'.join(map(str, new_fields)))

if __name__ == '__main__': 
  add_column()
