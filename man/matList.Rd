% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/matList.R
\name{matList}
\alias{matList}
\title{Generate a list of normalized matrices from BigWig files over genomic regions.}
\usage{
matList(
  bwf,
  bwr,
  names,
  grl,
  wins = list(Gene = 10),
  mode = "coverage",
  output = "norm.matrix",
  strand = "rev",
  smooth = FALSE,
  extend = 0,
  w,
  include_target = TRUE,
  target_ratio = 0.5,
  k = 10,
  keep
)
}
\arguments{
\item{bwf}{A list of \code{rtracklayer::BigWigFile} objects for forward strand data (e.g., output of \code{\link{importBWlist}}).}

\item{bwr}{A list of \code{rtracklayer::BigWigFile} objects for reverse strand data.}

\item{names}{A character vector specifying the names to assign to the resulting list of matrices.}

\item{grl}{A \code{GRangesList} where each element defines a set of genomic regions (e.g., genes, TSSs).}

\item{wins}{A named list specifying the number of bins (windows) for each region type in \code{grl}.
If multiple entries are provided, features are assumed to be combined.}

\item{mode}{Character string specifying the normalization mode for \code{normalizeToMatrix}, such as "coverage" or "mean".}

\item{output}{Output format: either \code{"norm.matrix"} (default) to return \code{normalizedMatrix} objects,
or \code{"matrix"} to return base R matrices.}

\item{strand}{Strand specificity: \code{"rev"} (default) for reverse alignment, \code{"for"} for forward, or \code{"no"} to ignore strand.}

\item{smooth}{Logical; if \code{TRUE}, applies smoothing in \code{normalizeToMatrix}.}

\item{extend}{Integer; number of base pairs to extend around each region.}

\item{w}{Optional integer; smoothing window size (used when \code{smooth = TRUE}).}

\item{include_target}{Logical; whether to include the target region in the matrix.}

\item{target_ratio}{Numeric between 0 and 1; controls the fraction of the matrix allocated to the target region.}

\item{k}{Integer; number of bins for the matrix (used only if \code{wins} has length 1).}

\item{keep}{(Unused) Optional parameter placeholder for future development.}
}
\value{
A named list of matrices, either as \code{EnrichedHeatmap::normalizedMatrix} objects
or base R numeric matrices, depending on the \code{output} parameter.
}
\description{
This function takes forward and reverse strand BigWig files and a list of genomic region sets,
then returns a list of normalized matrices representing signal intensities over those regions.
It supports flexible control over binning, strand handling, normalization mode, smoothing,
and how the matrix represents target regions.
}
\examples{
\dontrun{
# Load BigWig files and genomic regions
forward_bw <- importBWlist(c("sample1.f.bw", "sample2.f.bw"), c("sample1", "sample2"))
reverse_bw <- importBWlist(c("sample1.r.bw", "sample2.r.bw"), c("sample1", "sample2"))

gene_grl <- GRangesList(
  gene1 = GRanges("chr1", IRanges(1000, 2000), strand = "+", name = "gene1_a"),
  gene2 = GRanges("chr1", IRanges(3000, 4000), strand = "-", name = "gene2_b")
)

# Create normalized matrix list with 20 bins
matList(
  bwf = forward_bw,
  bwr = reverse_bw,
  names = c("sample1", "sample2"),
  grl = gene_grl,
  wins = list("Gene" = 20)
)
}
}
