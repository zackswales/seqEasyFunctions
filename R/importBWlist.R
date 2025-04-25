#' Import a list of BigWig files.
#'
#' This function imports multiple BigWig files into a list of \code{rtracklayer::BigWigFile} objects.
#'
#' @param bwf A character vector specifying the paths to the BigWig files.
#' @param names A character vector specifying the names to assign to the elements of the resulting list.
#' @param selection An optional \code{GRanges} or \code{RangesList} object specifying the genomic regions to import from each BigWig file. If \code{NULL} (default), the entire file is imported.
#'
#' @return A named list of \code{rtracklayer::BigWigFile} objects.
#'
#' @importFrom rtracklayer import BigWigFile BigWigSelection
#' @importFrom purrr future_map
#'
#' @examples
#' # Assuming you have BigWig files named "file1.bw" and "file2.bw"
#' # and a GRanges object named 'regions'
#' \dontrun{
#' bigwig_files <- c("file1.bw", "file2.bw")
#' file_names <- c("sample1", "sample2")
#'
#' # Import the entire files
#' bw_list_full <- importBWlist(bigwig_files, file_names)
#'
#' # Import only the specified regions
#' bw_list_selected <- importBWlist(bigwig_files, file_names, selection = regions)
#' }
#' @export
importBWlist <- function(bwf,names,selection=NULL){
  if(is.null(selection)){
    bwl <- bwf |> future_map(~import(.x, format = "BigWig"))
  }
  else{
    bwl <- bwf |> future_map(~import(.x, format = "BigWig", selection = BigWigSelection(selection)))
  }
  names(bwl) <- names
  bwl
}
