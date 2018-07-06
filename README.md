# CovFuzze Plot
A python script to plot gene coverage and IP peaks with standard deviations calculated from multiple replicates. 

usage: gene_peak_plot.py [-h] [-g GENE] [-o OUT] --bams BAMS [BAMS ...] --bed 
                         BED [--gtf GTF] [-p PEAKS] -l LABELS [LABELS ...]
                         [-n NSUBPLOTS] [--normalize]

## dependencies:
- numpy, tested with 1.11.1
- pysam, tested with 0.9.1
- pandas, tested with 0.18.1
- seaborn, tested with 0.7.1 
- matplotlib, tested with 2.0.2
- pybedtools, tested with 0.7.10 - requires bedtools available in path

## required arguments:
```
  -o OUT, --out OUT     output prefix (incl. dir) 
  --bams BAMS [BAMS ...] 
                        bam files (separated by spaces) 
  --bed BED             bed file for region(s) of interest  
  -l LABELS [LABELS ...], --labels LABELS [LABELS ...] 
                        labels associated with bams - if replicates, use same 
                        labels
```
## optional arguments:
```
  -h, --help            show this help message and exit 
  -g GENE, --gene GENE  gene name (default = GeneDoe) 
  --gtf GTF             gtf file for gene to plot exons as shaded regions 
  -p PEAKS, --peaks PEAKS 
                        bed file with peaks
  -n NSUBPLOTS, --nsubplots NSUBPLOTS 
                        number of subplots- bams will be split evenly based on 
                        the order given (default = 1) 
  --normalize           normalize by gene length/summed coverage (default = 
                        False)
```

## Example:
```
covfuzze.py -o $outdir/$prefix --bams 
    alignments/sample_1_Input.sample.star.bam alignments/sample_1_IP.sample.star.bam 
    alignments/sample_2_Input.sample.star.bam alignments/sample_2_IP.sample.star.bam 
    alignments/sample_3_Input.sample.star.bam alignments/sample_3_IP.sample.star.bam 
    --bed ${gene}_exons.bed 
    -l sample_Input sample_IP sample_Input sample_IP sample_Input sample_IP 
    --gene ${gene} -n 1 --normalize
```
