
Compute sample saturation curve by subsampling the dataset. 

Example toy dataset: 100 cells, 90% mapping rate

```bash
python saturation_curve.py -b test/alignments.bam -n 100 -r 0.9 -o output.tsv
python scripts/plot_curve.py  output.tsv saturation.png 
```
`output.tsv` - contains the sequencing saturation statistics for to 10 (`-s`) subsampling steps   
`saturation.png` - contains the plot of the saturation curve

![Saturation curve](saturation.png)

