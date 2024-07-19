#' Produce a standardised table using the Iterative Proportional Fitting method
#'
#' This ancillary function compute the standardised version of the input cross-tab using the Iterative Proportional Fitting method
#'
#' @param table Input cross-tabulation.
#' @param delta The desired level of accuracy
#' @param marginal.type The desired marginal type
#' @param custom.row.totals Custom row totals
#' @param custom.col.totals Custom col totals
#'
#' @return List containing the cross-tabulation in standardised format and number of iteration.
#'
#' @keywords internal
#'
#'
#'
standardize_table <- function(table, delta = 0.001, marginal.type = "average",
                              custom.row.totals = NULL, custom.col.totals = NULL) {
  # Function to determine desired totals based on type
  get.desired.totals <- function(table, dim, type, custom.totals) {
    if (!is.null(custom.totals)) {
      if (length(custom.totals) != dim(table)[dim]) {
        stop(paste("Length of custom totals does not match the number of", if(dim == 1) "rows" else "columns"))
      }
      return(custom.totals)
    }

    totals <- if (dim == 1) rowSums(table) else colSums(table)
    n <- length(totals)

    switch(type,
           "average" = rep(mean(totals), n),
           "percent" = rep(100 / n, n),
           stop("Invalid marginal type. Use 'average', 'percent', or provide custom totals.")
    )
  }

  # Get desired row and column totals
  desired.row.totals <- get.desired.totals(table, 1, marginal.type, custom.row.totals)
  desired.col.totals <- get.desired.totals(table, 2, marginal.type, custom.col.totals)

  # Initialization
  current.table <- table
  max.iterations <- 10000  # To prevent infinite loops
  iteration.count <- 0
  converged <- FALSE

  # Begin the iterative standardization procedure
  while (!converged && iteration.count < max.iterations) {
    previous.table <- current.table

    # Row standardization
    row.totals <- rowSums(current.table)
    for (i in 1:nrow(current.table)) {
      current.table[i, ] <- current.table[i, ] / row.totals[i] * desired.row.totals[i]
    }

    # Column standardization
    col.totals <- colSums(current.table)
    for (j in 1:ncol(current.table)) {
      current.table[, j] <- current.table[, j] / col.totals[j] * desired.col.totals[j]
    }

    # Check for convergence
    if (max(abs(current.table - previous.table)) < delta) {
      converged <- TRUE
    }
    iteration.count <- iteration.count + 1
  }

  if (!converged) {
    warning("Standardization did not converge within the maximum number of iterations.")
  }

  return(list(
    "table.stand" = current.table,
    "n.iterations" = iteration.count
  ))
}
