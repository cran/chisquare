#' Compute the standard error of Goodman-Kruskal tau for the input table
#'
#' This ancillary function compute the standard error of Goodman-Kruskal tau. In the code, some annotations
#' refer to equations from Bishop, Y. M., Fienberg, S. E., & Holland, P. W. (2007). Discrete Multivariate Analysis:
#' Theory and Practice. Springer. ISBN 9780387728056
#'
#' @param cont_table Input cross-tabulation.
#'
#' @return List with two components: the se of tau considering the row variable as dependent, and the se of tau
#' when the column variable is considered dependent
#'
#' @keywords internal
#'
#'
#'
# Function to compute standard errors for both directions of Goodman-Kruskal tau
compute_gk_tau_se <- function(cont_table) {
  # Ensure input is in matrix format for consistent calculation
  if(is.data.frame(cont_table)) {
    cont_table <- as.matrix(cont_table)
  }

  # Helper function to compute standard error for a single direction
  compute_single_se <- function(P) {
    # Calculate basic quantities needed for all steps
    N <- sum(P)
    P <- P / N
    pi_plus <- rowSums(P)
    p_plus_j <- colSums(P)

    # Calculate v using second expression from equation 11.3-22
    v <- sum(sapply(1:nrow(P), function(i) {
      sapply(1:ncol(P), function(j) {
        if(p_plus_j[j] > 0) (P[i,j] - pi_plus[i]*p_plus_j[j])^2/p_plus_j[j] else 0
      })
    }))

    # Calculate delta (denominator term)
    delta <- 1 - sum(pi_plus^2)

    # Initialize variance sum for equation 11.3-24
    var_sum <- 0

    # Calculate variance components
    for(i in 1:nrow(P)) {
      for(j in 1:ncol(P)) {
        if(P[i,j] == 0) next
        term1 <- 2 * v * sum(pi_plus[-i])
        term2 <- 2 * sum(P[-i,j]/p_plus_j[j])
        term3 <- sum(sapply(1:ncol(P), function(l) {
          if(p_plus_j[l] > 0) sum(P[-i,l]/p_plus_j[l]) else 0
        }))
        bracket_term <- term1 - delta * (term2 - term3)
        var_sum <- var_sum + (P[i,j] * bracket_term^2)
      }
    }

    # Return standard error with numerical stability check
    sqrt(max(var_sum / (N * delta^4), 0))
  }

  # Calculate and return standard errors for both directions
  # Column-dependent case uses original table
  # Row-dependent case uses transposed table

  return(
    list(
    row_dependent = compute_single_se(cont_table),
    col_dependent = compute_single_se(t(cont_table))
  ))
}
