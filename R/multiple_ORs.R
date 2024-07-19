#' Compute multiple ORs for tables larger than 2x2
#'
#' This ancillary function compute multiple ORs for tables larger than 2x2
#'
#' @param freq_table Input cross-tabulation.
#'
#' @return Table of multiple ORs
#'
#' @keywords internal
#'
#'
#'
# Function to calculate odds ratios using input table categories' labels or default labels if not provided;
# The odds ratios are computed for adjacent rows and columns, so representing
# the minimal set needed to describe all pairwise relationships in the contingency table.
calculate_odds_ratios <- function(freq_table) {
  # Check if the input is a matrix or data frame
  if(!is.matrix(freq_table) && !is.data.frame(freq_table)) {
    stop("Input must be a matrix or a data frame")
  }
  # Convert data frame to matrix if needed
  if(is.data.frame(freq_table)) {
    freq_table <- as.matrix(freq_table)
  }
  # Number of rows and columns
  n_rows <- nrow(freq_table)
  n_cols <- ncol(freq_table)
  # Get row and column names, if they are NULL, generate default labels
  row_labels <- rownames(freq_table)
  col_labels <- colnames(freq_table)
  if (is.null(row_labels)) {
    row_labels <- paste0("r", 1:n_rows)
  }
  if (is.null(col_labels)) {
    col_labels <- paste0("c", 1:n_cols)
  }
  # Initialize the odds ratios table
  odds_ratios <- matrix(NA, nrow = n_rows - 1, ncol = n_cols - 1)
  # Use row and column labels for the odds ratios table
  rownames(odds_ratios) <- paste(row_labels[1:(n_rows - 1)], "-", row_labels[2:n_rows], sep="")
  colnames(odds_ratios) <- paste(col_labels[1:(n_cols - 1)], "-", col_labels[2:n_cols], sep="")
  # Calculate the odds ratios
  for (i in 1:(n_rows - 1)) {
    for (j in 1:(n_cols - 1)) {
      # Extract the 2x2 subtable
      subtable <- freq_table[i:(i+1), j:(j+1)]

      # Check if either diagonal has a zero
      if (subtable[1,1] * subtable[2,2] == 0 || subtable[1,2] * subtable[2,1] == 0) {
        # Add 0.5 to all cells in the 2x2 subtable
        subtable <- subtable + 0.5
      }

      # Calculate the odds ratio using the (potentially) adjusted subtable
      odds_ratio <- (subtable[1,1] * subtable[2,2]) / (subtable[1,2] * subtable[2,1])
      odds_ratios[i, j] <- odds_ratio
    }
  }
  return(odds_ratios)
}
