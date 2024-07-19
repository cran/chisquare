#' Compute phi_max for the input table
#'
#' This ancillary function compute the maximum value ohi can achieve given the observed marginals
#'
#' @param crosstab Input cross-tabulation.
#'
#' @return Numerical value of phi_max.
#'
#' @keywords internal
#'
#'
#'
# Function to perform the required calculations based on a 2x2 cross-tabulation
compute_phi_max <- function(crosstab) {
  # Check if input is a dataframe and convert to matrix
  if (is.data.frame(crosstab)) {
    crosstab <- as.matrix(crosstab)
  }

  # Ensure the input is a 2x2 table
  if (!all(dim(crosstab) == c(2, 2))) {
    stop("Input must be a 2x2 table.")
  }

  # Calculate marginal totals
  row_totals <- margin.table(crosstab, 1) # Sum over columns (rows' totals)
  col_totals <- margin.table(crosstab, 2) # Sum over rows (columns' totals)

  # Find the maximum marginal total and identify the corresponding marginal
  if (max(row_totals) > max(col_totals)) {
    # Row marginal has the largest count
    largest_marginal <- row_totals
    other_marginal <- col_totals
  } else {
    # Column marginal has the largest count
    largest_marginal <- col_totals
    other_marginal <- row_totals
  }

  # Sort counts to find the smaller and larger counts within the identified marginal
  sorted_largest_marginal <- sort(largest_marginal)
  smaller_count_largest <- sorted_largest_marginal[1]
  larger_count_largest <- sorted_largest_marginal[2]

  # Perform the division within the identified marginal set
  division_result_largest <- smaller_count_largest / larger_count_largest

  # Sort other_counts to find the smaller and larger counts in the other marginal set
  sorted_other_marginal <- sort(other_marginal, decreasing = TRUE)
  larger_count_other <- sorted_other_marginal[1]
  smaller_count_other <- sorted_other_marginal[2]

  # Calculate the ratio of the larger count to the smaller count in the other marginal set
  division_result_other <- larger_count_other / smaller_count_other

  # Calculate the square root of the product of the two ratios
  final_result <- sqrt(division_result_largest * division_result_other)

  return(final_result)
}
