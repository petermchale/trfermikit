## Motivation 

Most SVs missed by short-read callers lie in tandem repeats: 

![](images/most_missing_SVs_lie_in_tandem_repeats.png)

(based upon Supplementary Data 53 of [Chaisson et al 2019](https://pubmed.ncbi.nlm.nih.gov/30992455) and private communication with Mark Chaisson) 

`trfermikit` has been optimized to finds SVs (currently only DELs) in a user-supplied set of tandem repeats. 

## `trfermikit` is more sensitive than `manta`

We assessed the performance of `trfermikit` 
and `manta`, in both cases relative to a long-read benchmark callset. Results can be found [here](evaluate-calls/evaluate.ipynb).

## How does it work?

`trfermikit` is based upon the [fermikit](https://pubmed.ncbi.nlm.nih.gov/26220959/) pipeline for deep Illumina resequencing data, which assembles reads into unitigs, maps them to the reference genome and then calls variants from the alignment.

`trfermikit` biases the minimap2 alignment step of the fermikit pipeline towards revealing deletions
by:

* increasing the reward for single-base matches
* increasing the penalty for single-base mismatches 
* decreasing the gap-open penalties (there are two because the cost function of gap length is piecewise linear)
* decreasing the gap-extension penalties 

This, by itself, recovers a lot of events that a more stringent caller would throw out. 
The trade-off is a large false-discovery rate. 

`trfermikit` mitigates this effect by throwing out calls that:
* are supported by "dirty" fermikit unitigs (essentially, those that have lots of small blocks when aligned to the reference or those whose mapping quality is zero)
* occur in “clusters”


## How fast is it?

If one confines the search to tandem repeats larger than 100bp in length:

```
trfermikit --min-repeat-length 100 [...]
```

then the pipeline takes < 1.5hr for a 70X genome.

## How does one use it?

An example of how to use `trfermikit` can be found [here](test-trfermikit.sh). 

```
conda create --name trfermikit python=3.8
pip install --requirements requirements.txt 
```

## TODO

Create a docker container and nextflow workflow and register both at [dockstore](https://dockstore.org/).

