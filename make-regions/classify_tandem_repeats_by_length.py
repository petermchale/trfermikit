import sys
import color_traceback
import argparse
import gzip 

def classify(fields, args): 
  _, start, end, period = fields
  condition1 = int(end) - int(start) > args.min_repeat_length
  condition2 = int(period) > args.min_repeat_period
  return 1 if condition1 and condition2 else 0

def add_column():
  parser = argparse.ArgumentParser(description='')
  parser.add_argument('--min-repeat-length', type=int, dest='min_repeat_length', help='')
  parser.add_argument('--min-repeat-period', type=int, dest='min_repeat_period', help='')  
  parser.add_argument('--log-file', type=str, dest='log_file', help='')  
  args = parser.parse_args()

  with gzip.open(args.log_file, 'wt') as log_file: 
    for line_number, line in enumerate(sys.stdin):
      if line_number == 0:
        if line != 'chrom\tchromStart\tchromEnd\tperiod\n': 
          from color_text import error
          error('repeats file does not have correct column names')
          sys.exit(1)
        print('{}\t{}'.format(line.strip(), 'class'), file=log_file)
      else: 
        fields = line.strip().split()
        classification = classify(fields, args)
        stdout_fields = fields[:3] + [classification]
        print('\t'.join(map(str, stdout_fields)))        
        log_fields = fields[:] + [classification]
        print('\t'.join(map(str, log_fields)), file=log_file)        

if __name__ == '__main__': 
  add_column()
