#' Generate a metaplot from a list of matrices.
#'
#' This function takes a list of matrices (e.g., coverage values across genomic regions)
#' and generates a metaplot, optionally comparing it to a control list of matrices.
#' It provides flexibility in grouping, faceting, and displaying error bars.
#'
#' @param matl A list of matrices, where each matrix represents a sample and rows
#'   correspond to genomic positions.
#' @param matlc An optional list of control matrices with the same dimensions and
#'   structure as `matl`. If provided, the function will plot the log2 ratio of
#'   `matl` to `matlc`. Defaults to `NULL`.
#' @param feature A string specifying the label for the x-axis (e.g., "Gene",
#'   "Transcript"). Defaults to "Gene".
#' @param unit A string specifying the label for the y-axis (e.g., "Coverage (BPM)",
#'   "Signal Intensity"). Defaults to "Coverage (BPM)".
#' @param title A string specifying the title of the plot. Defaults to "Gene metaplot".
#' @param breaks A numeric vector specifying the positions of the x-axis breaks.
#'   Defaults to `c(0, 20, 60, 80)`.
#' @param labels A character vector specifying the labels for the x-axis breaks.
#'   The length of `labels` must be the same as the length of `breaks`.
#'   Defaults to `c("-200", "TSS", "TTS", "+200")`.
#' @param colmap A named vector of colors to use for different groups (defined by
#'   `colour_by`). The names of the vector should correspond to the unique values
#'   in the column specified by `colour_by`.
#' @param split An optional data frame whose rownames correspond to the rownames
#'   of the matrices in `matl`. This data frame can contain columns to group the
#'   data before plotting (e.g., by experimental condition). Defaults to `NULL`.
#' @param facet An optional string or vector of strings specifying the column(s)
#'   in `split` to use for faceting the plot. Defaults to `NULL`.
#' @param angle A numeric value specifying the angle (in degrees) for the x-axis
#'   labels. Defaults to `0`.
#' @param strip_fill A string specifying the fill color for the facet strip
#'   background. Defaults to "white".
#' @param facet_scale A string specifying whether the scales should be fixed
#'   ("fixed"), free ("free"), free in the x direction ("free_x"), or free in
#'   the y direction ("free_y"). Passed to `facet_wrap` or `facet_grid2`.
#'   Defaults to "fixed".
#' @param max_quantile A numeric value between 0 and 1 specifying the upper
#'   quantile to use for outlier removal when calculating the mean or sum.
#'   Defaults to `1` (no upper outlier removal).
#' @param min_quantile A numeric value between 0 and 1 specifying the lower
#'   quantile to use for outlier removal when calculating the mean or sum.
#'   Defaults to `0` (no lower outlier removal).
#' @param pseudo A numeric value added to both `matl` and `matlc` before log2
#'   transformation to avoid taking the logarithm of zero. Only used when
#'   `matlc` is provided. Defaults to `1`.
#' @param alpha A numeric value between 0 and 1 specifying the alpha (transparency)
#'   of the lines. Defaults to `1` (fully opaque).
#' @param linewidth A numeric value specifying the width of the lines. Defaults
#'   to `0.5`.
#' @param error A logical value indicating whether to plot error bars (standard
#'   error of the mean). Only applicable when `summarise_by = "mean"`.
#'   Defaults to `FALSE`.
#' @param alpha_error A numeric value between 0 and 1 specifying the alpha
#'   (transparency) of the error bar ribbons. Defaults to `0.5`.
#' @param facet_nrow An integer specifying the number of rows to use when
#'   `facet_type = "wrap"`. Defaults to `2`.
#' @param facet_type A string specifying the faceting type. Can be "wrap"
#'   (using `facet_wrap`) or "grid" (using `facet_grid2`). Defaults to "wrap".
#' @param facet_independent A logical value indicating whether the scales should
#'   be independent for each panel in `facet_grid2`. Only used when
#'   `facet_type = "grid"`. Defaults to `FALSE`.
#' @param summarise_by A string specifying how to summarise the data within
#'   groups. Can be "mean" or "sum". Defaults to "mean".
#' @param colour_by A string specifying the column name in the combined data frame
#'   (or in `split` if provided) to use for coloring the lines. Defaults to
#'   "Sample".
#'
#' @return A ggplot object representing the metaplot.
#'
#' @import ggplot2
#' @import dplyr
#' @import purrr
#' @import tidyr
#' @import stringr
#' @importFrom rlang .data sym !! !!! :=
#' @importFrom matrixStats rowSds
#' @importFrom ggplot2 facet_wrap facet_grid2
#'
#' @examples
#' # Create some dummy data
#' set.seed(123)
#' mat1 <- matrix(rnorm(100), nrow = 10)
#' mat2 <- matrix(rnorm(100, 0.5), nrow = 10)
#' matl_example <- list(sample1 = mat1, sample2 = mat2)
#' rownames(mat1) <- rownames(mat2) <- paste0("gene", 1:10)
#'
#' # Basic metaplot
#' mplot(matl_example)
#'
#' # Metaplot with custom title and axis labels
#' mplot(matl_example, title = "My Metaplot", feature = "Region", unit = "Signal")
#'
#' # Metaplot with custom colors
#' col_map <- c("sample1" = "blue", "sample2" = "red")
#' mplot(matl_example, colmap = col_map)
#'
#' # Add a splitting factor
#' split_df <- data.frame(condition = rep(c("A", "B"), each = 5),
#'                        row.names = paste0("gene", 1:10))
#' mplot(matl_example, split = split_df, colour_by = "condition")
#'
#' # Facet the plot
#' mplot(matl_example, split = split_df, facet = "condition")
#'
#' # Metaplot with control data
#' matlc_example <- list(sample1 = matrix(rnorm(100, 0.1), nrow = 10),
#'                       sample2 = matrix(rnorm(100, 0.2), nrow = 10))
#' rownames(matlc_example[[1]]) <- rownames(matlc_example[[2]]) <- paste0("gene", 1:10)
#' mplot(matl_example, matlc = matlc_example, unit = "log2(Treatment/Control)")
#'
#' # Metaplot with error bars
#' mplot(matl_example, error = TRUE)
#'
#' @export
mplot<-function(matl,matlc=NULL,feature = "Gene", unit = "Coverage (BPM)", title = "Gene metaplot", breaks = c(0,20,60,80), labels = c("-200","TSS","TTS","+200"), colmap, split=NULL, facet=NULL, angle = 0, strip_fill = "white",facet_scale="fixed", max_quantile = 1, min_quantile = 0, pseudo = 1,alpha=1, linewidth=0.5, error = F, alpha_error=0.5,facet_nrow=2,facet_type="wrap",facet_independent=F,summarise_by="mean",colour_by="Sample"){

  ml <- matl
  if(!is.null(matlc)){
    ml <- map2(matl,matlc,~log2((.x + pseudo) / (.y + pseudo)))
  }

  df <- ml |> imap(~as.data.frame(.x) |>
                     mutate(Sample = .y) |>
                     rownames_to_column("name")) |>
    bind_rows()

  if(!is.null(split)){
    if(error == T & summarise_by == "mean"){
      df.error <- df |>
        left_join(split |> rownames_to_column("name"),by="name") |>
        group_by(Sample,across(all_of(names(split)))) |>
        summarise(across(-name,~std.error(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T))) |>
        ungroup()
    }
    if(summarise_by == "mean"){
      df <- df |>
        left_join(split |> rownames_to_column("name"),by="name") |>
        group_by(Sample,across(all_of(names(split)))) |>
        summarise(across(-name,~mean(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T))) |>
        ungroup()
    }
    if(summarise_by == "sum"){
      df <- df |>
        left_join(split |> rownames_to_column("name"),by="name") |>
        group_by(Sample,across(all_of(names(split)))) |>
        summarise(across(-name,~sum(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T))) |>
        ungroup()
    }
  }
  else{
    if(error == T & summarise_by == "mean"){
      df.error <- df |> group_by(Sample) |>
        summarise(across(-name,~std.error(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T)))
    }
    df <- df |> group_by(Sample) |>
      summarise(across(-name,~mean(.x[.x <= quantile(.x,max_quantile) & .x >= quantile(.x,min_quantile)],na.rm=T)))
  }

  names <- c("Sample")
  if(!is.null(split)){
    names <- c(names,names(split))
  }

  names(df) <- c(names,paste0("w",1:(ncol(df) - length(names))))
  df <- df |> pivot_longer(-all_of(names),names_to = "Index",values_to = "Coverage") |>
    mutate(Index = Index |> str_remove("w") |> as.numeric())

  if(error == T & summarise_by == "mean"){
    names(df.error) <- c(names,paste0("w",1:(ncol(df.error) - length(names))))
    df.error <- df.error |> pivot_longer(-all_of(names),names_to = "Index",values_to = "Error") |>
      mutate(Index = Index |> str_remove("w") |> as.numeric())
    df <- df |> left_join(df.error,by=join_by(!!!names,Index))
  }


  p <- df |> ggplot(aes(Index,Coverage,group=!!sym(colour_by),colour=!!sym(colour_by)))

  if(error == T & summarise_by == "mean"){
    p <- p +geom_ribbon(aes(ymin=Coverage-Error,ymax=Coverage+Error,fill=!!sym(colour_by)),alpha=alpha_error)
  }

  p <- p +  #geom_vline(xintercept = vlines, colour = "darkgrey") +
    geom_line(alpha=alpha,linewidth=linewidth) +
    theme_bw() +
    labs(x = feature,y = unit, title = title) +
    scale_colour_manual(values = colmap) +
    scale_fill_manual(values = colmap) +
    scale_x_continuous(breaks=breaks,labels=labels,minor_breaks = NULL) +
    theme(axis.text.x = element_text(angle = angle, hjust=1),
          strip.background = element_rect(fill=strip_fill))

  if(!is.null(facet)){
    if(facet_type=="wrap"){
      p <- p + facet_wrap(facet,scales = facet_scale,nrow = facet_nrow)
    }
    else{
      p <- p + facet_grid2(facet,scales = facet_scale,independent=facet_independent)
    }
  }

  p
}
