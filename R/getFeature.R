#' Get genomic ranges based on specified features and flanking regions.
#'
#' This function takes a GRanges object and extracts or creates new genomic
#' ranges based on specified start and end features (TSS, TES, or Exon) and
#' their flanking regions. It allows for specifying exon numbers and boundaries
#' for precise range definition.
#'
#' @param object A \code{GRanges} object containing genomic features (e.g., genes
#'        with exons). This object must have metadata columns including 'blockStarts',
#'        'blockSizes', and 'blockCount' if 'Exon' is specified as a feature.
#' @param start_feature A character string specifying the starting feature.
#'        Must be one of "TSS" (Transcription Start Site), "TES" (Transcription
#'        End Site), or "Exon". Default is "TSS".
#' @param start_flank An integer indicating the number of bases to flank the
#'        \code{start_feature}. Default is 0.
#' @param start_exon An integer specifying the exon number to use if
#'        \code{start_feature} is "Exon".
#' @param start_exon_boundary A character string specifying the exon boundary to
#'        use if \code{start_feature} is "Exon". Must be either "5prime" or
#'        "3prime".
#' @param start_direction A character string indicating the direction of the
#'        flank relative to the \code{start_feature}. Must be either "up"
#'        (upstream) or "down" (downstream). Default is "up".
#' @param end_feature A character string specifying the ending feature.
#'        Must be one of "TSS", "TES", or "Exon". Default is "TES".
#' @param end_flank An integer indicating the number of bases to flank the
#'        \code{end_feature}. Default is 0.
#' @param end_direction A character string indicating the direction of the
#'        flank relative to the \code{end_feature}. Must be either "up"
#'        (upstream) or "down" (downstream). Default is "down".
#' @param end_exon An integer specifying the exon number to use if
#'        \code{end_feature} is "Exon".
#' @param end_exon_boundary A character string specifying the exon boundary to
#'        use if \code{end_feature} is "Exon". Must be either "5prime" or
#'        "3prime".
#' @param return A character string specifying which ranges to return. Must be
#'        one of "all" (both plus and minus strand ranges), "plus" (only plus
#'        strand ranges), or "minus" (only minus strand ranges). Default is "all".
#'
#' @return A \code{GRanges} object containing the extracted or created genomic ranges.
#'
#' @import GenomicRanges
#'
#' @examples
#' # Create a dummy GRanges object
#' library(GenomicRanges)
#' genes <- GRanges(
#'   seqnames = "chr1",
#'   ranges = IRanges(start = c(100, 2000), end = c(1000, 3000)),
#'   strand = c("+", "-"),
#'   blockStarts = c("0,200,500", "0,300"),
#'   blockSizes = c("100,300,500", "200,500"),
#'   blockCount = c(3, 2)
#' )
#'
#' # Get TSS regions with 100bp upstream flank
#' tss_regions <- getFeature(genes, start_feature = "TSS", start_flank = 100)
#' tss_regions
#'
#' # Get regions from TSS to TES
#' gene_bodies <- getFeature(genes, start_feature = "TSS", end_feature = "TES")
#' gene_bodies
#'
#' # Get regions around the 2nd exon start
#' exon2_start <- getFeature(genes, start_feature = "Exon", start_exon = 2,
#'                            start_exon_boundary = "5prime", start_flank = 50)
#' exon2_start
#'
#' # Get regions spanning the 1st to the 2nd exon
#' exon1_to_2 <- getFeature(genes, start_feature = "Exon", start_exon = 1,
#'                           start_exon_boundary = "5prime",
#'                           end_feature = "Exon", end_exon = 2,
#'                           end_exon_boundary = "3prime")
#' exon1_to_2
#'
#' @export
getFeature <- function(object, start_feature = "TSS", start_flank = 0, start_exon = NULL,
                       start_exon_boundary = NULL, start_direction = "up", end_feature = "TES",
                       end_flank = 0, end_direction = "down", end_exon = NULL, end_exon_boundary = NULL,
                       return = "all")
{
  if (!inherits(object, "GRanges")) {
    stop("Input must be a GRanges object.")
  }
  allowedFeatures = c("TSS", "TES", "Exon")
  if (!start_feature %in% allowedFeatures) {
    stop("start_feature must be either 'TSS', 'TES' or 'Exon'.")
  }
  if (!end_feature %in% allowedFeatures) {
    stop("end_feature must be either 'TSS', 'TES' or 'Exon'.")
  }
  allowed_directions = c("up", "down")
  if (!start_direction %in% allowed_directions) {
    stop("start_direction should be set to either 'up' or 'down'.")
  }
  if (!end_direction %in% allowed_directions) {
    stop("stop_direction should be set to either 'up' or 'down'")
  }
  starts = list()
  object = object
  plusStrandRanges = GRanges()
  minusStrandRanges = GRanges()
  if (start_feature == "Exon") {
    if (!is.null(start_exon_boundary)) {
      exonStartNumber = start_exon
      splitIndex = object$blockCount >= exonStartNumber
      if (TRUE %in% splitIndex) {
        object = object[splitIndex]
        if (end_feature != "Exon") {
          plusStrandRanges <- subset(object, strand(object) ==
                                       "+")
          minusStrandRanges <- subset(object, strand(object) ==
                                        "-")
        }
        else {
          if (!is.null(end_exon)) {
            exonEndNumber = end_exon
            splitIndex = object$blockCount >= exonEndNumber
            if (TRUE %in% splitIndex) {
              object = object[splitIndex]
              plusStrandRanges <- subset(object, strand(object) ==
                                           "+")
              minusStrandRanges <- subset(object, strand(object) ==
                                            "-")
            }
            else {
              stop("There are no genes with that many exons. Please enter another exon end.")
            }
          }
          else {
            stop("If end_feature is set to Exon then the end_exon must not be null.")
          }
        }
      }
      else {
        stop("There are no genes with that many exons. Please enter another exon start.")
      }
    }
    else {
      stop("If the start feature is set to Exon, start_exon_boundary must be set to either '3prime' or '5prime'.")
    }
  }
  else if (start_feature == "TSS") {
    if (is.null(start_exon) & is.null(start_exon_boundary)) {
      if (end_feature != "Exon") {
        plusStrandRanges <- subset(object, strand(object) ==
                                     "+")
        minusStrandRanges <- subset(object, strand(object) ==
                                      "-")
      }
      else {
        if (!is.null(end_exon)) {
          exonEndNumber = end_exon
          splitIndex = object$blockCount >= exonEndNumber
          if (TRUE %in% splitIndex) {
            object = object[splitIndex]
            plusStrandRanges <- subset(object, strand(object) ==
                                         "+")
            minusStrandRanges <- subset(object, strand(object) ==
                                          "-")
          }
          else {
            stop("There are no genes with that many exons. Please enter another exon end.")
          }
        }
        else {
          stop("If end_feature is set to Exon then end_exon must not be null")
        }
      }
    }
    else {
      stop("If start feature is TSS then no exon start information must be entered.")
    }
  }
  else {
    if (is.null(start_exon) & is.null(start_exon_boundary)) {
      if (end_feature != "Exon") {
        plusStrandRanges <- subset(object, strand(object) ==
                                     "+")
        minusStrandRanges <- subset(object, strand(object) ==
                                      "-")
      }
      else {
        if (!is.null(end_exon)) {
          exonEndNumber = end_exon
          splitIndex = object$blockCount >= exonEndNumber
          if (TRUE %in% splitIndex) {
            object = object[splitIndex]
            plusStrandRanges <- subset(object, strand(object) ==
                                         "+")
            minusStrandRanges <- subset(object, strand(object) ==
                                          "-")
          }
          else {
            "There are no genes with that many exons. Please enter another exon end."
          }
        }
        else {
          stop("If end_feature is set to Exon then end_exon must not be null")
        }
      }
    }
    else {
      stop("If start feature is TES then no exon start information must be entered.")
    }
  }
  minusStrandTempRanges = minusStrandRanges |> as.data.frame()
  plusStrandTempRanges = plusStrandRanges |> as.data.frame()
  if (start_feature == "TSS") {
    if (start_direction == "up") {
      plusStrandTempRanges$start = start(plusStrandRanges) -
        start_flank
      minusStrandTempRanges$end = end(minusStrandRanges) +
        start_flank
    }
    else {
      plusStrandTempRanges$start = start(plusStrandRanges) +
        start_flank

      minusStrandTempRanges$end = end(minusStrandRanges) -
        start_flank
    }
  }
  if (start_feature == "TES") {
    if (start_direction == "up") {
      plusStrandTempRanges$start = end(plusStrandRanges) -
        start_flank
      minusStrandTempRanges$end = start(minusStrandRanges) +
        start_flank
    }
    else {
      plusStrandTempRanges$start = end(plusStrandRanges) +
        start_flank
      minusStrandTempRanges$end = start(minusStrandRanges) -
        start_flank
    }
  }
  if (start_feature == "Exon") {
    if (start_exon_boundary == "5prime") {
      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts,
                                           ","), as.numeric)
      exon_starts_plus <- vector(mode = "numeric", length = length(plusStrandRanges))
      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] = tss_plus[i] + block_starts_plus[[i]][start_exon]
      }
      start_minus = start(minusStrandRanges)
      block_counts_minus = minusStrandRanges$blockCount
      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts,
                                           ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus,
      )
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes,
                                          ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus,
      )
      exon_starts_minus = vector(mode = "numeric", length = length(minusStrandRanges))
      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][start_exon] +
          block_sizes_minus[[i]][start_exon]
      }
    }
    else {
      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount
      block_sizes_plus = lapply(strsplit(plusStrandRanges$blockSizes,
                                         ","), as.numeric)
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts,
                                           ","), as.numeric)
      exon_starts_plus <- vector(mode = "numeric", length = length(plusStrandRanges))
      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] <- tss_plus[i] + block_starts_plus[[i]][start_exon] +
          block_sizes_plus[[i]][start_exon]
      }
      start_minus = start(minusStrandRanges)
      block_counts_minus = minusStrandRanges$blockCount
      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts,
                                           ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus,
      )
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes,
                                          ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus,
      )
      exon_starts_minus = vector(mode = "numeric", length = length(minusStrandRanges))
      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][start_exon]
      }
    }
    if (start_direction == "up") {
      plusStrandTempRanges$start = exon_starts_plus -
        start_flank
      minusStrandTempRanges$end = exon_starts_minus +
        start_flank
    }
    else {
      plusStrandTempRanges$start = exon_starts_plus +
        start_flank
      minusStrandTempRanges$end = exon_starts_minus -
        start_flank
    }
  }
  if (end_feature == "TSS") {
    if (end_direction == "up") {
      plusStrandTempRanges$end = start(plusStrandRanges) -
        end_flank
      minusStrandTempRanges$start = end(minusStrandRanges) +
        end_flank
    }
    else {
      plusStrandTempRanges$end = start(plusStrandRanges) +
        end_flank
      minusStrandTempRanges$start = end(minusStrandRanges) -
        end_flank
    }
  }
  if (end_feature == "TES") {
    if (end_direction == "up") {
      plusStrandTempRanges$end = end(plusStrandRanges) -
        end_flank
      minusStrandTempRanges$start = start(minusStrandRanges) +
        end_flank
    }
    else {
      plusStrandTempRanges$end = end(plusStrandRanges) +
        end_flank
      minusStrandTempRanges$start = start(minusStrandRanges) -
        end_flank
    }
  }
  else if (end_feature == "Exon") {
    if (end_exon_boundary == "5prime") {
      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts,
                                           ","), as.numeric)
      exon_starts_plus <- vector(mode = "numeric", length = length(plusStrandRanges))
      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] = tss_plus[i] + block_starts_plus[[i]][end_exon]
      }
      start_minus = start(minusStrandRanges)
      block_counts_minus = minusStrandRanges$blockCount
      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts,
                                           ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus,
      )
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes,
                                          ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus,
      )
      exon_starts_minus = vector(mode = "numeric", length = length(minusStrandRanges))
      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][end_exon] +
          block_sizes_minus[[i]][end_exon]
      }
    }
    else {
      tss_plus <- start(plusStrandRanges)
      block_counts_plus <- plusStrandRanges$blockCount
      block_sizes_plus = lapply(strsplit(plusStrandRanges$blockSizes,
                                         ","), as.numeric)
      block_starts_plus <- lapply(strsplit(plusStrandRanges$blockStarts,
                                           ","), as.numeric)
      exon_starts_plus <- vector(mode = "numeric", length = length(plusStrandRanges))
      for (i in seq_along(plusStrandRanges)) {
        exon_starts_plus[i] <- tss_plus[i] + block_starts_plus[[i]][end_exon] +
          block_sizes_plus[[i]][end_exon]
      }
      start_minus = start(minusStrandRanges)
      block_counts_minus = minusStrandRanges$blockCount
      block_starts_minus = lapply(strsplit(minusStrandRanges$blockStarts,
                                           ","), as.numeric)
      block_starts_minus = revElements(block_starts_minus,
      )
      block_sizes_minus = lapply(strsplit(minusStrandRanges$blockSizes,
                                          ","), as.numeric)
      block_sizes_minus = revElements(block_sizes_minus,
      )
      exon_starts_minus = vector(mode = "numeric", length = length(minusStrandRanges))
      for (i in seq_along(minusStrandRanges)) {
        exon_starts_minus[i] <- start_minus[i] + block_starts_minus[[i]][end_exon]
      }
    }
    if (end_direction == "up") {
      plusStrandTempRanges$end = exon_starts_plus - end_flank
      minusStrandTempRanges$start = exon_starts_minus +
        end_flank
    }
    else {
      plusStrandTempRanges$end = exon_starts_plus + end_flank
      minusStrandTempRanges$start = exon_starts_minus -
        end_flank
    }
  }
  plusStrandTempRanges <- GRanges(plusStrandTempRanges)
  minusStrandTempRanges <- GRanges(minusStrandTempRanges)
  if (return == "all") {
    allRanges = c(minusStrandTempRanges, plusStrandTempRanges)
    return(allRanges)
  }
  else if (return == "plus") {
    return(plusStrandTempRanges)
  }
  else if (return == "minus") {
    return(minusStrandTempRanges)
  }
  else {
    print("Please set a return value to 'all', 'plus' or 'minus'")
  }
}
