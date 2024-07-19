#' Compute the chi-squared-maximising version of the input table
#'
#' This ancillary function computes the chi-squared-maximising version of the input table
#'
#' @param observed_table Input cross-tabulation.
#'
#' @return List containing the chi-squared-maximising version of the input table and the chi-squared statistic
#'
#' @keywords internal
#'
#'
#'
maximize_chi_squared <- function(observed_table) {
  # Store original row and column names
  row_names <- rownames(observed_table)
  col_names <- colnames(observed_table)

  # Ensure the input is a matrix
  observed_table <- as.matrix(observed_table)

  row_totals <- rowSums(observed_table)
  col_totals <- colSums(observed_table)
  max_table <- matrix(0, nrow = nrow(observed_table), ncol = ncol(observed_table))

  while (any(row_totals > 0) && any(col_totals > 0)) {
    for (i in seq_along(row_totals)) {
      for (j in seq_along(col_totals)) {
        if (row_totals[i] == col_totals[j] && row_totals[i] > 0) {
          max_table[i, j] <- row_totals[i]
          row_totals[i] <- 0
          col_totals[j] <- 0
        }
      }
    }
    largest_row_total <- which.max(row_totals)
    largest_col_total <- which.max(col_totals)
    if (row_totals[largest_row_total] > 0 && col_totals[largest_col_total] > 0) {
      smaller_of_two <- min(row_totals[largest_row_total], col_totals[largest_col_total])
      max_table[largest_row_total, largest_col_total] <- smaller_of_two
      row_totals[largest_row_total] <- row_totals[largest_row_total] - smaller_of_two
      col_totals[largest_col_total] <- col_totals[largest_col_total] - smaller_of_two
    }
  }

  expected <- outer(rowSums(observed_table), colSums(observed_table), FUN = "*") / sum(observed_table)
  chi_squared_max <- sum((max_table - expected)^2 / expected)

  # Restore row and column names to max_table
  rownames(max_table) <- row_names
  colnames(max_table) <- col_names

  return(list("max_table" = max_table, "chi_squared_max" = chi_squared_max))
}

