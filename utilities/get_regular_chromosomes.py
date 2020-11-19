import sys 
import color_traceback

chromosomes = set('chr{}'.format(id) for id in list(range(1,23)) + ['X', 'Y'])

for line in sys.stdin: 
  chromosome = line.strip().split()[0]
  if chromosome in chromosomes: 
    print(line, end="") 

