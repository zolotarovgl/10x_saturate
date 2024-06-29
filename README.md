
Compute sample saturation curve by subsampling the dataset. 

Example toy dataset: 100 cells, 90% mapping rate

```bash
python saturation_curve.py -b test/alignments.bam -n 100 -r 0.9 -o output.tsv 
```
`output.tsv` contains the sequencing saturation statistics for to 10 (`-s`) subsampling steps   

## TODOS:  

[] - saturation prediction script  


