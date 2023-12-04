#' Suggest chi-squared testing method for the input contingency table
#'
#' This function returns a suggested method for chi-square testing, on the basis of different criteria related to the
#' features of the input cross-tabulation.
#'
#' @param cross_tab Input cross-tabulation.
#'
#' @return A carachter vector containing a suggestion for the chi-square testing method to employ.
#'
#' @keywords internal
#'
#' @importFrom stats pchisq chisq.test
#'
#'
suggest_chi_squared_method <- function(cross_tab) {
  # Store the original setting
  original_warn <- options("warn")

  # Deactivate warning messages
  options(warn = -1)

  # Calculate the grand total and expected counts
  grand_total <- sum(cross_tab)
  num_cells <- nrow(cross_tab) * ncol(cross_tab)
  expected <- round(chisq.test(cross_tab)$expected, 3)
  min_expected <- min(expected)

  # Check if the cross-tab is 2x2 or larger
  cross_tab_size <- ifelse(nrow(cross_tab) == 2 && ncol(cross_tab) == 2, "2x2", "larger")

  # Decision logic with explanations
  explanation <- ""
  if (cross_tab_size == "2x2") {
    explanation <- paste0(explanation, "Since the table is 2x2, ")
    if (grand_total >= 5 * num_cells) {
      explanation <- paste0(explanation, "and the grand total (", grand_total,") is equal to or larger than 5 times the number of cells (", 5 * num_cells, "), use the traditional Chi-squared test. Permutation or Monte Carlo method can also be considered")
    } else {
      if (min_expected >= 1) {
        explanation <- paste0(explanation, "the grand total (", grand_total,") is smaller than 5 times the number of cells (", 5 * num_cells, "), and the minimum expected count (", min_expected, ") is equal to or larger than 1, use the (N-1)/N adjusted Chi-squared test. Permutation or Monte Carlo method can also be considered")
      } else {
        explanation <- paste0(explanation, "the grand total (", grand_total,") is smaller than 5 times the number of cells (", 5 * num_cells, "), and the minimum expected count (", min_expected, ") is less than 1, use the Permutation or Monte Carlo method")
      }
    }
  } else { # Larger than 2x2
    explanation <- paste0(explanation, "Since the table is larger than 2x2, ")
    if (grand_total >= 5 * num_cells) {
      explanation <- paste0(explanation, "and the grand total (", grand_total,") is equal to or larger than 5 times the number of cells (", 5 * num_cells, "), use the traditional Chi-squared test. Permutation or Monte Carlo method can also be considered")
    } else {
      if (min_expected >= 1) {
        explanation <- paste0(explanation, "the grand total (", grand_total,") is smaller than 5 times the number of cells (", 5 * num_cells, "), and the minimum expected count (", min_expected, ") is equal to or larger than 1, use the (N-1)/N adjusted Chi-squared test. Permutation or Monte Carlo method can also be considered")
      } else {
        explanation <- paste0(explanation, "the grand total (", grand_total,") is smaller than 5 times the number of cells (", 5 * num_cells, "), and the minimum expected count (", min_expected, ") is less than 1, use the Permutation or Monte Carlo method")
      }
    }
  }

  # Restore the original warning setting
  options(warn = original_warn$warn)

  return(explanation)
}
