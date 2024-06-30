
## Compute sample saturation curve by subsampling the dataset. 

This script subsamples the .bam file and fits the saturation curve based on the amount of __input reads__. For the estimation  to be precise, the mapping rate and the number of cells should be supplied.  

Example toy dataset: 1000 cells, 70% mapping rate

```bash
# downsample the reads 
python saturation_table.py -b test/sample.bam -n 1000 -r 0.7 -o output.tsv
# fit the MM model, predict the number of input reads and plot
python scripts/plot_curve.py  output.tsv saturation.png 
```
`output.tsv` - contains the sequencing saturation statistics for to 10 (`-s`) subsampling steps   
`saturation.png` - contains the plot of the saturation curve

![Saturation curve](img/saturation.png)

# Speed-up   

The saturation can in principle be estimated from the reads coming from a subset of chromosomes.  

__TODO: add speed vs estimated saturation plot__  
 

