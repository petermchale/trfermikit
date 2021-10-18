import json 

counts = {}
for region_size_range in [
  "600,625", "625,650", "650,700", "700,800", "800,1000", "1000,2000", "2000,100000"
]:
  region_counts_filename = f'HGSVC2/data/GWD,HG02818/truvari-DEL-{region_size_range}-pacbio-manta/summary.txt'
  with open(region_counts_filename, 'r') as fh: 
    region_counts = json.load(fh)
    counts[region_size_range] = {
      'TP': region_counts['TP-base'],
      'FP': region_counts['FP'],
      'FN': region_counts['FN']
    }
with open('aaron.json', 'w') as fh:     
  json.dump(counts, fh, indent=2)


