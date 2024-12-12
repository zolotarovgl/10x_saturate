
## Compute sample saturation curve by downsampling 

__TODOs:__  
- check if supplying the .bam file with unaligned reads makes a difference  
- handle different .bam files - with mapped and unmapped reads - better reporting     

This script subsamples the .bam file and fits the saturation curve based on the amount of __input reads__. 
For the estimation  to be precise, the mapping rate and the number of cells should be supplied.  
Example toy dataset: 1000 cells, 80% mapping rate
Downsample the reads, extract the tags and compute saturation stats (use `-c` to use more cores): 

```bash
python saturation_table.py --ncpu 5 -b test/sample.bam -n 1000 -r 0.8 -o output.tsv
```
__CAVE: this creates a big tabular file__ 

Fit the MM model, predict the number of input reads for `--target` saturation and plot:
```bash
python scripts/plot_curve.py  output.tsv saturation.png --target 0.7 
```
`output.tsv` - contains the sequencing saturation statistics for to 10 (`-s`) subsampling steps   
`saturation.png` - contains the plot of the saturation curve

![Saturation curve](img/saturation.png)

It's useful to examine residuals plot to see if the model tends to over or underestimate the coverage needeed.   

![Residuals](img/saturation_residuals.png)

## How to compute saturation?  

The [saturation](https://kb.10xgenomics.com/hc/en-us/articles/115005062366-What-is-sequencing-saturation) is defined as following:  
> Sequencing saturation is a measure of the fraction of library complexity that was sequenced in a given experiment. The inverse of one minus the sequencing saturation can be interpreted as the number of additional reads it would take to detect a new transcript.


How is it calculated?  


Vis [10x reference](https://kb.10xgenomics.com/hc/en-us/articles/115003646912-How-is-sequencing-saturation-calculated) for details.

> The web_summary.html output from cellranger count includes a metric called "Sequencing Saturation". This metric quantifies the fraction of reads originating from an already-observed UMI. More specifically, this is the fraction of confidently mapped, valid cell-barcode, valid UMI reads that are non-unique (match an existing cell-barcode, UMI, gene combination).

The script attampts to mimic the `cellranger`'s approach by extracting the relevant tags from the `.bam` file. Note that `cellranger` is computing saturation by considering confidently mapped reads which may inroduce some discrepancies with other methods. 


The formula for calculating this metric is as follows:

$$ \text{Sequencing Saturation} = 1 - \frac{N_{\text{dedup}}}{N_{\text{dedup}} + N_{\text{dup}}}$$

