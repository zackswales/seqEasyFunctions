#' Generate a list of enriched heatmaps from a list of matrices.
#'
#' This function takes a list of matrices (typically the output of \code{\link{matList}})
#' and generates a list of enriched heatmaps using the \code{EnrichedHeatmap} package.
#' It allows for customization of color schemes, row splitting, column features,
#' quantile-based scaling, row name display, window labels, y-axis limits,
#' summarisation of enrichment, axis labels, row clustering, and log2 transformation.
#'
#' @param matl A named list of matrices, where each matrix represents signal intensity
#'        over genomic regions (e.g., output of \code{\link{matList}}).
#' @param wins A named integer vector specifying the number of bins (windows) that the
#'        matrices in \code{matl} were normalized to for each feature type. The names
#'        should correspond to the feature types (e.g., \code{c("Gene" = 10)}).
#' @param split An optional \code{data.frame} used to split the rows of the heatmaps.
#'        Each column in the data frame will define a way to split the rows.
#' @param split_cols A named list specifying the colors for the levels of the splitting variable
#'        provided in \code{split}. The names of the list should correspond to the column
#'        names in the data frame used for splitting.
#' @param max_quantile A numeric value between 0 and 1 specifying the upper quantile for
#'        scaling the color map. Values above this quantile will be capped.
#' @param min_quantile A numeric value between 0 and 1 specifying the lower quantile for
#'        scaling the color map. Values below this quantile will be capped.
#' @param col_fun A character string specifying the color function to use. Options include
#'        "red" (white to red), "bl2rd" (blue to white to red), or "red0" (white to light red to dark red starting from 0).
#'        Alternatively, a custom color function generated by \code{circlize::colorRamp2} can be provided.
#' @param show_row_names A logical value indicating whether to display row names on the heatmaps.
#' @param win_labels An optional character vector specifying the labels for the feature windows
#'        displayed in the top annotation. If \code{NULL}, the names of \code{wins} are used.
#' @param ylim An optional numeric vector of length 2 specifying the y-axis limits for the
#'        enrichment profile in the top annotation. If \code{NULL}, limits are automatically determined.
#' @param summarise_by A character string specifying the function to summarise the rows for the
#'        enrichment profile in the top annotation (e.g., "mean", "median").
#' @param axis_labels A character string specifying the label for the x-axis of the
#'        enrichment profile in the top annotation.
#' @param row_km An integer specifying the number of clusters to perform k-means clustering on the rows.
#'        If \code{NULL} or 1, no clustering is performed.
#' @param log2 A logical value indicating whether to apply a log2 transformation (with a pseudocount of 1)
#'        to the matrices in \code{matl} before plotting.
#'
#' @return A list of \code{EnrichedHeatmap} objects, one for each matrix in the input \code{matl}.
#'
#' @importFrom EnrichedHeatmap EnrichedHeatmap HeatmapAnnotation anno_enriched
#' @importFrom circlize colorRamp2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices gpar
#' @importFrom stats quantile
#' @importFrom purrr imap
#' @importFrom graphics unit
#'
#' @examples
#' # Assuming you have generated a list of matrices using matList()
#' \dontrun{
#' # Create a dummy matrix list for demonstration
#' mat_list <- list(
#'   sample1 = matrix(rnorm(100), nrow = 20),
#'   sample2 = matrix(rnorm(100), nrow = 20)
#' )
#' wins_info <- c("FeatureA" = 5)
#'
#' # Generate a list of enriched heatmaps
#' hm_list <- hmList(matl = mat_list, wins = wins_info)
#'
#' # Generate heatmaps with row splitting and custom colors
#' split_df <- data.frame(group = factor(rep(c("Group1", "Group2"), each = 10)))
#' split_colors <- list("group" = c("Group1" = "blue", "Group2" = "green"))
#' hm_list_split <- hmList(matl = mat_list, wins = wins_info, split = split_df,
#'                         split_cols = split_colors)
#' }
#' @export
hmList <- function(matl, wins, split = NULL, split_cols, max_quantile = 0.99, min_quantile = 0, col_fun = "red", show_row_names = TRUE, win_labels = NULL, ylim = NULL, summarise_by = "mean", axis_labels = "", row_km = 1, log2 = FALSE) {

  ## ColourMap
  reds <- RColorBrewer::brewer.pal(n = 9, name = "Reds")

  if (log2) {
    matl <- lapply(matl, function(x) log2(x + 1))
  }

  common_min <- quantile(unlist(matl), min_quantile)
  common_max <- quantile(unlist(matl), max_quantile)
  if (col_fun == "red") {
    col_fun <- circlize::colorRamp2(c(common_min, common_max), c("white", "red"))
  } else if (col_fun == "bl2rd") {
    col_fun <- circlize::colorRamp2(c(common_min, 0, common_max), c("blue", "white", "red"))
  } else if (col_fun == "red0") {
    col_fun <- circlize::colorRamp2(c(0, 1, common_max), c("white", reds[3], reds[7]))
  }

  if (is.null(ylim)) {
    ymin <- min(0, common_min)
    ymax <- common_max
    offset <- quantile(1:(ymax - ymin), 0.01)
    ylim <- c(ymin - offset, ymax + offset)
  }

  features <- vector()
  fcols <- RColorBrewer::brewer.pal(9, "Set1")[1:length(wins)]
  names(fcols) <- names(wins)
  for (i in 1:length(wins)) {
    features <- append(features, rep(names(wins)[i], wins[i]))
  }

  if (is.null(win_labels)) {
    win_labels <- names(wins)
  }

  if (!is.null(split)) {

    rowAnno <- rowAnnotation(df = split, col = split_cols, show_legend = FALSE)

    hml <- matl |> imap(~EnrichedHeatmap(
      .x,
      name = .y,
      column_title = .y,
      col = col_fun,
      show_row_names = show_row_names,
      row_names_gp = gpar(fontsize = 5),
      axis_name = axis_labels,
      row_split = split,
      row_km = row_km,
      column_title_gp = gpar(fontsize = 10),
      cluster_rows = FALSE,
      top_annotation = c(
        HeatmapAnnotation(
          Features = features,
          show_legend = TRUE,
          border = TRUE,
          col = list(Features = fcols),
          show_annotation_name = FALSE,
          annotation_legend_param = list(at = names(wins), labels = win_labels)
        ),
        HeatmapAnnotation(
          enriched = anno_enriched(
            value = summarise_by,
            gp = gpar(fontsize = 4, lwd = 2, col = split_cols[[1]]),
            axis_param = list(side = "left", facing = "inside"),
            ylim = ylim
          )
        ),
        gap = unit(2, "mm")
      )
    ))
    hml$rowAnno <- rowAnno
    hml
  } else {
    hml <- matl |> imap(~EnrichedHeatmap(
      .x,
      name = .y,
      column_title = .y,
      col = col_fun,
      show_row_names = show_row_names,
      row_names_gp = gpar(fontsize = 5),
      axis_name = axis_labels,
      column_title_gp = gpar(fontsize = 10),
      row_km = row_km,
      cluster_rows = FALSE,
      top_annotation = c(
        HeatmapAnnotation(
          Features = features,
          show_legend = TRUE,
          border = TRUE,
          col = list(Features = fcols),
          show_annotation_name = FALSE,
          annotation_legend_param = list(at = names(wins), labels = win_labels)
        ),
        HeatmapAnnotation(
          enriched = anno_enriched(
            gp = gpar(fontsize = 4, lwd = 2),
            axis_param = list(side = "left", facing = "inside"),
            ylim = ylim
          )
        ),
        gap = unit(2, "mm")
      )
    ))
    hml
  }
}
