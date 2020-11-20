import sys 
import color_traceback

def classify(fields):
  tandem_repeat_labels = fields[-1] 
  return max(map(int, tandem_repeat_labels.split(','))) 

def add_column(): 
  for line in sys.stdin: 
    old_fields = line.strip().split()
    new_fields = old_fields + [classify(old_fields)]
    print('\t'.join(map(str, new_fields)))

if __name__ == '__main__': 
  add_column()
