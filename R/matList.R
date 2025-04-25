#' Generate a list of normalized matrices from BigWig files over genomic regions.
#'
#' This function takes forward and reverse strand BigWig files and a list of genomic region sets,
#' then returns a list of normalized matrices representing signal intensities over those regions.
#' It supports flexible control over binning, strand handling, normalization mode, smoothing,
#' and how the matrix represents target regions.
#'
#' @param bwf A list of \code{rtracklayer::BigWigFile} objects for forward strand data (e.g., output of \code{\link{importBWlist}}).
#' @param bwr A list of \code{rtracklayer::BigWigFile} objects for reverse strand data.
#' @param names A character vector specifying the names to assign to the resulting list of matrices.
#' @param grl A \code{GRangesList} where each element defines a set of genomic regions (e.g., genes, TSSs).
#' @param wins A named list specifying the number of bins (windows) for each region type in \code{grl}.
#' If multiple entries are provided, features are assumed to be combined.
#' @param mode Character string specifying the normalization mode for \code{normalizeToMatrix}, such as "coverage" or "mean".
#' @param output Output format: either \code{"norm.matrix"} (default) to return \code{normalizedMatrix} objects,
#' or \code{"matrix"} to return base R matrices.
#' @param strand Strand specificity: \code{"rev"} (default) for reverse alignment, \code{"for"} for forward, or \code{"no"} to ignore strand.
#' @param smooth Logical; if \code{TRUE}, applies smoothing in \code{normalizeToMatrix}.
#' @param extend Integer; number of base pairs to extend around each region.
#' @param w Optional integer; smoothing window size (used when \code{smooth = TRUE}).
#' @param include_target Logical; whether to include the target region in the matrix.
#' @param target_ratio Numeric between 0 and 1; controls the fraction of the matrix allocated to the target region.
#' @param k Integer; number of bins for the matrix (used only if \code{wins} has length 1).
#' @param keep (Unused) Optional parameter placeholder for future development.
#'
#' @return A named list of matrices, either as \code{EnrichedHeatmap::normalizedMatrix} objects
#' or base R numeric matrices, depending on the \code{output} parameter.
#'
#' @importFrom rtracklayer import
#' @importFrom GenomicRanges strand subset
#' @importFrom purrr future_map map imap bind_cols
#' @importFrom MatrixGenerics rowMeans2
#' @importFrom EnrichedHeatmap normalizeToMatrix as.normalizedMatrix
#'
#' @examples
#' \dontrun{
#' # Load BigWig files and genomic regions
#' forward_bw <- importBWlist(c("sample1.f.bw", "sample2.f.bw"), c("sample1", "sample2"))
#' reverse_bw <- importBWlist(c("sample1.r.bw", "sample2.r.bw"), c("sample1", "sample2"))
#'
#' gene_grl <- GRangesList(
#'   gene1 = GRanges("chr1", IRanges(1000, 2000), strand = "+", name = "gene1_a"),
#'   gene2 = GRanges("chr1", IRanges(3000, 4000), strand = "-", name = "gene2_b")
#' )
#'
#' # Create normalized matrix list with 20 bins
#' matList(
#'   bwf = forward_bw,
#'   bwr = reverse_bw,
#'   names = c("sample1", "sample2"),
#'   grl = gene_grl,
#'   wins = list("Gene" = 20)
#' )
#' }
#' @export
matList <- function(bwf, bwr, names, grl, wins = list("Gene" = 10), mode = "coverage", output = "norm.matrix",
                    strand = "rev", smooth = FALSE, extend = 0, w, include_target = TRUE, target_ratio = 0.5,
                    k = 10, keep) {

  fbw <- bwf
  #fbw <- bwf |> future_map(~import(.x, format = "BigWig"))
  #names(fbw) <- names

  if(strand %in% c("for", "rev")) {
    rbw <- bwr
    #rbw <- bwr |> future_map(~import(.x, format = "BigWig"))
    #names(rbw) <- names

    fgrl = grl |> map(~subset(.x, strand(.x) == "+"))
    rgrl = grl |> map(~subset(.x, strand(.x) == "-"))
  }

  if(length(wins) > 1) { ## Length of features to combine
    matl <- names |> future_map(function(s) {

      if(strand == "rev") {
        fmat <- fgrl |> imap(~normalizeToMatrix(rbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y],
                                                mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |>
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)
        rmat <- rgrl |> imap(~normalizeToMatrix(fbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y],
                                                mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |>
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)

        mat <- rbind(fmat, rmat)
        rownames(mat) <- c(fgrl[[1]]$name, rgrl[[1]]$name)
      }
      else if(strand == "for") {
        fmat <- fgrl |> imap(~normalizeToMatrix(fbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y],
                                                mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |>
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)
        rmat <- rgrl |> imap(~normalizeToMatrix(rbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y],
                                                mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |>
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)

        mat <- rbind(fmat, rmat)
        rownames(mat) <- c(fgrl[[1]]$name, rgrl[[1]]$name)
      }
      else if(strand == "no") {
        mat <- grl |> imap(~normalizeToMatrix(fbw[[s]], .x, extend = 0, value_column = "score", k = wins[.y],
                                              mean_mode = mode, smooth = smooth)) |> bind_cols() |> as.data.frame() |>
          as.matrix() |> as.normalizedMatrix(k_target = sum(wins), extend = 0)
        rownames(mat) <- grl[[1]]$name
      }

      mat

    }, .progress = TRUE)
  }
  else {  ## Single feature - regular normalizeToMatrix
    matl <- names |> future_map(function(s) {

      if(strand == "rev") {
        fmat <- normalizeToMatrix(rbw[[s]], fgrl[[1]], extend = extend, value_column = "score", k = wins[1],
                                  mean_mode = mode, smooth = smooth, w = w, include_target = include_target,
                                  target_ratio = target_ratio)
        rmat <- normalizeToMatrix(fbw[[s]], rgrl[[1]], extend = extend, value_column = "score", k = wins[1],
                                  mean_mode = mode, smooth = smooth, w = w, include_target = include_target,
                                  target_ratio = target_ratio)

        mat <- rbind(fmat, rmat)
        rownames(mat) <- c(fgrl[[1]]$name, rgrl[[1]]$name)
      }
      else if(strand == "for") {
        fmat <- normalizeToMatrix(fbw[[s]], fgrl[[1]], extend = extend, value_column = "score", k = wins[.y],
                                  mean_mode = mode, smooth = smooth, w = w, include_target = include_target,
                                  target_ratio = target_ratio)
        rmat <- normalizeToMatrix(rbw[[s]], rgrl[[1]], extend = extend, value_column = "score", k = wins[1],
                                  mean_mode = mode, smooth = smooth, w = w, include_target = include_target,
                                  target_ratio = target_ratio)

        mat <- rbind(fmat, rmat)
        rownames(mat) <- c(fgrl[[1]]$name, rgrl[[1]]$name)
      }
      else if(strand == "no") {
        mat <- normalizeToMatrix(fbw[[s]], grl[[1]], extend = extend, value_column = "score", k = wins[.y],
                                 mean_mode = mode, smooth = smooth, w = w, include_target = include_target,
                                 target_ratio = target_ratio)
        rownames(mat) <- grl[[1]]$name
      }

      mat

    }, .progress = TRUE)
  }

  names(matl) <- names

  if(output == "norm.matrix") {
    matl
  }
  else if(output == "matrix") {
    matl |> map(~as.data.frame(.x) |> as.matrix())
  }
}
