#' Ancillary function for median polish analysis
#'
#' This ancillary function compute the standardised Pearson median polish residual and their adjusted version,
#' following Simonoff 2003 and Mosteller-Parunak 1985.
#'
#' @param df Input cross-tabulation.
#'
#' @return List containing robust expected counts, standardised and adjusted standardised Pearson median polish residuals
#'
#' @keywords internal
#'
#'
#'
compute_median_polish_residuals <- function(df) {
  # Add a small constant to avoid log(0)
  log_counts <- log(as.matrix(df) + 0.001)

  # Apply median polish
  median_polish_result <- medpolish(log_counts, trace.iter = FALSE)

  # Calculate robust expected values
  log_expected <- outer(median_polish_result$row,
                        median_polish_result$col, "+") +
    median_polish_result$overall
  robust_expected <- exp(log_expected)

  # Calculate Pearson median polish residuals
  pearson_mp_residuals <- (as.matrix(df) - robust_expected) /
    sqrt(robust_expected)

  # Calculate Haberman-adjusted median polish residuals
  n <- sum(df)
  row_props <- rowSums(df)/n
  col_props <- colSums(df)/n

  # Calculate adjustment factors for each cell
  adj_factors <- matrix(0, nrow = nrow(df), ncol = ncol(df))
  for (i in 1:nrow(df)) {
    for (j in 1:ncol(df)) {
      adj_factors[i,j] <- sqrt((1 - row_props[i]) * (1 - col_props[j]))
    }
  }

  # Apply adjustment to residuals
  haberman_mp_residuals <- pearson_mp_residuals / adj_factors

  # Set dimension names to match the input data
  dimnames(pearson_mp_residuals) <- dimnames(df)
  dimnames(haberman_mp_residuals) <- dimnames(df)
  dimnames(robust_expected) <- dimnames(df)

  # Return all the calculated values
  return(list(
    "robust_expected" = robust_expected,
    "pearson_mp_residuals" = pearson_mp_residuals,
    "haberman_mp_residuals" = haberman_mp_residuals
  ))
}
