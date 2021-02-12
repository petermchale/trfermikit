## Motivation 

Most SVs missed by short-read callers lie in tandem repeats: 

![](images/most_missing_SVs_lie_in_tandem_repeats.png)

(based upon Supplementary Data 53 of [Chaisson et al 2019](https://pubmed.ncbi.nlm.nih.gov/30992455) and private communication with Mark Chaisson). 
This observation prompted us to attempt to discover novel SVs in tandem repeats. 

## Scope 

As callers exist 
to capture SVs in tandem repeats where 
the repeat unit is smaller than 6bps,
known in the community as Short Tandem Repeats 
(STRs),
we designed trfermikit to pick up SVs that manta missed
in tandem repeats with repeat units 
larger than 6bps, known as 
Variable Number Tandem Repeats (VNTRs). 

## Impact 

We assessed the performance of trfermikit and manta, in both cases relative to a long-read benchmark callset, on VNTRs. 
The results show that [(a)](experiments/paper_figures/TPRs_FDRs/DEL.svg) trfermikit has better sensitivity-FDR trade-offs than manta
and [(b)](experiments/paper_figures/manta_complementarity/DEL-manta.svg) trfermikit is complementary to manta. (In the figures, the red circle indicates the default operating point of trfermikit.)

## How does it work?

trfermikit is based upon the [fermikit](https://pubmed.ncbi.nlm.nih.gov/26220959/) pipeline for deep Illumina resequencing data, which assembles reads into unitigs, maps them to the reference genome, and then calls variants from the alignment.

trfermikit biases the `minimap2` alignment step of the fermikit pipeline towards revealing deletions
by:
* increasing the reward for single-base matches
* increasing the penalty for single-base mismatches 
* decreasing the gap-open penalties (there are two because the cost function of gap length is piecewise linear)
* decreasing the gap-extension penalties 

This recovers a lot of events that a more stringent caller would throw out at the cost of an elevated number of false discoveries. 
trfermikit mitigates this by throwing out calls that:
* are supported by "dirty" fermikit unitigs (essentially, those that have lots of small blocks when aligned to the reference or those whose mapping quality is zero)
* occur in “clusters”

## Installation

```
git clone https://github.com/petermchale/trfermikit
cd trfermikit
bash install.sh 
conda activate trfermikit
```
Only installation on Linux x86_64 is currently supported.

## Usage 

Assuming that the path to this directory on your filesystem is 
`${root}`, usage is: 

```
PATH="${root}:$PATH"

trfermikit [OPTIONS] \
  --output <path to results directory> \
  --reference <path to reference fasta (without the ".fa" suffix)> \
  --alignments <path to short-read alignments (without the ".cram" suffix"; index assumed to be present)> \
  --threads <number of threads>
```

Options are: 
``` 
--hg19 
      Use hg19 build of the human reference genome. 
      If this flag is not specified, trfermikit uses build hg38.
--functional-regions PATH 
      Restrict examination to those tandem repeats 
      that lie in the regions indicated by PATH
      (without the ".bed.gz" suffix). 
      [default value: None]
--min-repeat-length INT
      Only consider tandem repeats 
      whose total number of bps is larger than INT 
      [default value: 0].
```


## Output 

Discovered SVs are output in indexed vcf format to the results directory at
```
fermikit.raw.decomposed.normalized.${svtype}.unitigSupport.thinned.vcf.gz
```
where `${svtype}` is either `DEL` or `INS`. 

Tandem-repeat regions used to discover SVs appear in the results directory in indexed bed format at
```
regions.bed.gz
```

The parameter configuration used to make the discoveries appear in the results directory at `config.json`.



## How fast is it?

If one confines the search to tandem repeats with total lengths larger than 100bp 
(where one finds most tandem-repeat-associated DELs larger than 50bp),

```
trfermikit --min-repeat-length 100 ...
```

then the pipeline takes about 2 hours for a 70X genome.


## TODO

* Create a docker container and nextflow workflow (with nextflow processes for "make-regions", "make-calls" and "filter-calls") and register both at [dockstore](https://dockstore.org/).
* Make `--genome-build` and `--threads` optional, with default values of "hg38" and "1", respectively 
* Support bams

