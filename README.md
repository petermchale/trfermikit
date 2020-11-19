## Motivation 

Most SVs missed by short-read callers lie in tandem repeats: 

![](images/most_missing_SVs_lie_in_tandem_repeats.png)

(based upon Supplementary Data 53 of Chaisson et al 2019 and private communication with Mark Chaisson) 

`trfermikit` has been optimized to finds SVs (currently only DELs) in a user-supplied set of tandem repeats. 

## How it works 

`trfermikit` is based upon the [fermikit](https://github.com/lh3/fermikit) pipeline for deep Illumina resequencing data, which assembles reads into unitigs, maps them to the reference genome and then calls variants from the alignment.

`trfermikit` biases the alignment step of the fermikit pipeline towards revealing deletions:

```
minimap2 -A10 -B12 -O6,26 -E1,0 [...]
```

and optimizes the variant calling step as follows 
(-V relates to "per-base divergence"; -q relates to mapping quality):

```
htsbox pileup -V1 -q1
```

## Usage

An example of how to use `trfermikit` can be found [here](test-trfermikit.sh). 

## Runtime

If one confines the search to tandem repeats larger than 100bp in length:

```
trfermikit --min-repeat-length 100 [...]
```

then the pipeline takes < 1.5hr for a 70X genome

## Evaluation 

We assessed the performance of `trfermikit` 
and `manta` (in both cases relative to a long-read benchmark callset). Results can be found [here](evaluate-calls/evaluate.ipynb).

## TODO

Create a docker container and nextflow workflow and register both at [dockstore](https://dockstore.org/).

