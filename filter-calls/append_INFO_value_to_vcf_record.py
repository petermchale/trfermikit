import sys
from color_text import error 
import color_traceback

INFO_index = 7 

key = sys.argv[1]

for line in sys.stdin:
  fields = line.strip().split('\t')
  INFO = fields[INFO_index]
  if INFO == '.':
    error('INFO does not contain key: {}'.format(key))
    sys.exit(1)
  key_values = INFO.split(';') 
  INFO = dict([key_value.split('=') for key_value in key_values])
  print('\t'.join(fields) + '\t{}\n'.format(INFO[key]), end='')
