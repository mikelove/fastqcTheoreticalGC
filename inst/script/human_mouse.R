library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

set.seed(1)
head(seqlengths(Hsapiens),25)
generateDistn(genome=Hsapiens, nchrom=25,
              file="fastqc_theoretical_gc_hg38_genome.txt",
              name="Human Genome (UCSC hg38)")

set.seed(1)
hs.txs <- extractTranscriptSeqs(Hsapiens, TxDb.Hsapiens.UCSC.hg38.knownGene)
generateDistn(seqs=hs.txs,
              file="fastqc_theoretical_gc_hg38_txome.txt",
              name="Human Transcriptome (UCSC hg38)")

library(BSgenome.Mmusculus.UCSC.mm10)
library(TxDb.Mmusculus.UCSC.mm10.knownGene)

set.seed(1)
head(seqlengths(Mmusculus),22)
generateDistn(genome=Mmusculus, nchrom=22,
              file="fastqc_theoretical_gc_mm10_genome.txt",
              name="Mouse Genome (UCSC mm10)")

set.seed(1)
mm.txs <- extractTranscriptSeqs(Mmusculus, TxDb.Mmusculus.UCSC.mm10.knownGene)
generateDistn(seqs=mm.txs,
              file="fastqc_theoretical_gc_mm10_txome.txt",
              name="Mouse Transcriptome (UCSC mm10)")
