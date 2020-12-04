import sys 
import color_traceback

def classify(fields): 
  # TODO:
  # convert comment line to dictionary for each line, and filter there on period size

  min_length = sys.argv[1]
  _, start, end = fields
  return 1 if int(end) - int(start) > int(min_length) else 0

def add_column(): 
  for line in sys.stdin: 
    old_fields = line.strip().split()
    new_fields = old_fields + [classify(old_fields)]
    print('\t'.join(map(str, new_fields)))

if __name__ == '__main__': 
  add_column()
