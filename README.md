# fastqcTheoreticalGC

R package for generating theoretical GC content curves for the FASTQC module

The `generateDistn` function will create a theoretical distribution
from a genome or transcriptome, and then write this out to a file
which can be read by MultiQC, to be plotted over the GC content 
curves. The function randomly samples n=1 million simulated reads
(default 100bp) to generate the theoretical distribution.

See `?generateDistn` for more details on the function.

See the example file `inst/script/human_mouse.R` for examples
of generating human and mouse genome and transcriptome theoretical
distributions (the output files are also included there).

<https://github.com/mikelove/fastqcTheoreticalGC/blob/master/inst/script/human_mouse.R>

An example of generating a transcriptome theoretical distribution is:

```{r}
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
set.seed(1)
hs.txs <- extractTranscriptSeqs(Hsapiens, TxDb.Hsapiens.UCSC.hg38.knownGene)
generateDistn(seqs=hs.txs,
              file="fastqc_theoretical_gc_hg38_txome.txt",
              name="Human Transcriptome (UCSC hg38)")
```

# Assumptions

The function generates reads uniformly from the sequences or chromosomes 
provided. So, for the transcriptome, the example files involve generation of 
reads from all the transcripts (with counts proportional to the length 
of the transcripts). Genes with more isoforms are contributing more reads. 
You could pick a single isoform per gene, or you could use the weights to 
generate reads proportional to estimated transcript expression from an 
RNA-seq experiment.

The reads are generated from all positions of the sequences (allowing 
for the read length). They do not take into account fragment length 
distribution or any kind of positional bias.

For the transcriptome / RNA-seq, the read starts do not take into 
account random hexamer bias, but are uniform across the transcripts.

From my personal research into DNA-seq and RNA-seq sequence bias, 
I believe that these simplifying assumptions are not of practical 
consquence for coming up with a theoretical distribution of GC content 
for quality control testing purposes.
If however, you want to estimate a theoretical distribution for RNA-seq 
which takes all of these other factors into acccount, you can use the 
Biconductor package [alpine](http://bioconductor.org/packages/alpine),
which implements the methods of [this paper](http://www.nature.com/nbt/journal/v34/n12/full/nbt.3682.html),
and the `plotGC` function which can generate these distributions, conditional
on sample-specific transcript expression, positional bias, 
random hexamer priming bias, etc.
