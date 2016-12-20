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
