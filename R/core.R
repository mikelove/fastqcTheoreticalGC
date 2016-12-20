#' Generate theoretical GC content distributions
#'
#' This function generates random simulated reads from
#' either provided \code{seqs} (best for RNA-seq)
#' or from a genome (best for DNA-seq). The GC content
#' of these reads is then tabulated to produce a distribution
#' file which can be read by MultiQC to be displayed
#' on top of the FASTQC GC content module. Either \code{seqs}
#' or \code{genome} is required, and only one can be specified.
#' Specifying \code{genome} requires also specifying \code{nchrom}.
#' 
#' @param seqs a DNAStringSet of the sequences to simulate
#' read from. E.g. for RNA-seq, the transcripts, which
#' can be generated with \code{extractTranscriptSeqs}
#' from the GenomicFeatures package.
#' See the example script located in \code{inst/script/human_mouse.R}
#' @param genome a BSgenome object.
#' See the example script located in \code{inst/script/human_mouse.R}
#' @param nchrom the number of chromosomes from the genome to simulate
#' reads from
#' @param file the path of the file to write out
#' @param n the number of reads to simulate
#' @param bp the basepair of the reads
#' @param wts optional weights to go along with the \code{seqs} or
#' the chromosomes in \code{genome}, e.g. to represent
#' more realistic expression of transcripts
#' @param name the name to be printed at the top of the file
#'
#' @return the name of the file which was written
#'
#' @references
#'
#' MultiQC: 
#' http://multiqc.info/
#' 
#' FASTQC:
#' http://www.bioinformatics.babraham.ac.uk/projects/fastqc/
#' 
#' @export
generateDistn <- function(seqs, genome, nchrom,
                          file="fastqc_theoretical_gc.txt",
                          n=1e6, bp=100, wts=1, name="") {
  stopifnot(missing(seqs) | missing(genome))
  if (!missing(genome)) stopifnot(!missing(nchrom))

  # first the routine generating GC from given sequences (txome)
  if (!missing(seqs)) {
    # remove seqs that are too short
    keep <- width(seqs) >= bp
    if (length(wts) > 1) {
      stopifnot(length(wts) == length(seqs))
      stopifnot(sum(keep) > 0)
      wts <- wts[keep]
    }
    message(paste(sum(keep),"sequences of sufficient length"))
    seqs <- seqs[keep]
    # weights (optional) can specify e.g. transcript expression
    prob <- wts * width(seqs)
    prob <- prob/sum(prob)
    idx <- sample(length(seqs), n, replace=TRUE, prob=prob)
    message(paste("generating",n,"reads"))
    molecules <- seqs[idx]
    # sample start positions uniformly
    starts <- round(runif(length(molecules),1,width(molecules)-bp+1))
    reads <- subseq(molecules,start=starts,width=bp)

  } else {
    # now the routine generating from genome
    chrom.lens <- head(seqlengths(genome),nchrom)
    message(paste0("will generate reads from ",names(chrom.lens)[1],
                  "-",names(chrom.lens)[nchrom]))
    stopifnot(all(chrom.lens > bp))
    if (length(wts) > 1) stopifnot(length(wts) == nchrom)
    prob <- chrom.lens * wts
    prob <- prob/sum(prob)
    chroms <- sample(names(chrom.lens),n,replace=TRUE,prob=prob)
    starts <- round(runif(n,1,chrom.lens[chroms]-bp+1))
    gr <- GRanges(chroms,IRanges(starts,width=bp))
    message(paste("generating",n,"reads"))
    # sorting speeds up accessing the genome
    gr <- sort(gr)
    reads <- as(Views(genome, gr),"DNAStringSet")
  }
  
  # tabulate frequencies per percentile of GC content
  gc <- as.vector(letterFrequency(reads,"GC"))
  total <- as.vector(letterFrequency(reads,"ACGT"))
  gc.content <- (gc/total)[total > 0]
  message(paste0("mean (sd) GC content: ",100*round(mean(gc.content),2),
                  " (",100*round(sd(gc.content),2),")"))
  dens <- hist(round(gc.content,2), breaks=0:101/100 - 1/200, plot=FALSE)$density
  dens <- round(dens, 3)
  out <- cbind(GC=0:100, dens)
  message(paste("writing out density to",file))
  if (name != "") {
    name <- paste(":",name)
  }
  write(paste0("# FastQC theoretical GC content curve",name), file=file)
  write.table(out, file=file, append=TRUE, sep="\t",
              quote=FALSE, row.names=FALSE, col.names=FALSE)
  file
}
