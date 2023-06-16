#' Visualize Odds Ratios for a 2xk (where k >= 2) Contingency Table (internal function)
#'
#' This function creates a plot of odds ratios with 95% confidence intervals
#' for a 2xk (where k >= 2) contingency table.
#'
#' @param ctable A 2xk (where k >= 2) contingency table as a matrix or data frame with row and column names.
#' @param reference.level The index of the reference level for odds ratio calculations (default: 1). The user must select the column level to serve as the reference level.
#' @param row.level The index of the row category to be used in odds ratio calculations (1 or 2). The user must select the row level to which the calculation of the odds ratios make reference (default: 1).
#' @param or.alpha The significance level used for the confidence intervals (default: 0.05).
#'
#' @return A plot of odds ratios with 95% confidence intervals for a 2xk contingency table.
#'
#' @keywords internal
#'
#' @importFrom graphics axis lines text
#'
#'
visualize_odds_ratios <- function(ctable, reference.level = 1, row.level = 1, or.alpha = 0.05) {
  # Convert data frame to matrix if needed
  if (is.data.frame(ctable)) {
    ctable <- as.matrix(ctable)
  }

  if (!any(dim(ctable) == 2) || ncol(ctable) < 2) {
    stop("To plot the odds ratios, the input cross-tab must have dimension 2xk (where k >= 2).
         You can either enter a different cross-tab or subset the existing one.")
  }

  if (reference.level > ncol(ctable)) {
    stop(paste0("Invalid reference.level. It must be between 1 and ", ncol(ctable), "."))
  }

  if (row.level != 1 && row.level != 2) {
    stop("Invalid row.level. It must be either 1 or 2.")
  }

  # Extract row.level name
  rc_name <- rownames(ctable)[row.level]

  # Extract level names
  level_names <- colnames(ctable)

  if (is.null(level_names)) {
    level_names <- 1:ncol(ctable)
  }

  n_levels <- length(level_names)

  # Calculate odds ratios and confidence intervals
  odds_ratios <- c()
  lower_bounds <- c()
  upper_bounds <- c()

  for (i in 1:n_levels) {
    if (i != reference.level) {
      # Check for zeros along the diagonal of the 2x2 table
      if (ctable[1, i] * ctable[2, reference.level] == 0 || ctable[1, reference.level] * ctable[2, i] == 0) {
        # Add 0.5 to every cell of the 2x2 table
        ctable_sub <- ctable[c(row.level, 3 - row.level), c(i, reference.level)] + 0.5
      } else {
        ctable_sub <- ctable[c(row.level, 3 - row.level), c(i, reference.level)]
      }

      # Corrected odds_ratio calculation
      odds_ratio <- (ctable_sub[1, 1] * ctable_sub[2, 2]) / (ctable_sub[1, 2] * ctable_sub[2, 1])

      se <- sqrt(sum(1 / ctable_sub))
      z <- qnorm(1 - or.alpha / 2)

      lower_bounds <- append(lower_bounds, exp(log(odds_ratio) - z * se))
      upper_bounds <- append(upper_bounds, exp(log(odds_ratio) + z * se))
    } else {
      odds_ratio <- 1
      lower_bounds <- append(lower_bounds, NA)
      upper_bounds <- append(upper_bounds, NA)
    }
    odds_ratios <- append(odds_ratios, odds_ratio)
  }


  # Create the visualization
  plot(NULL, xlim = c(0.5, n_levels + 0.5), ylim = c(0, max(upper_bounds, na.rm = TRUE) * 1.2),
       xlab = "", ylab = "Odds Ratios", xaxt = 'n', yaxt = 'n',
       main = "Odds Ratios with 95% Confidence Intervals", cex.main = 0.85,
       sub = paste0("The ORs are relative to the '", rc_name, "' row level",
                    "\nReference column level: '", level_names[reference.level], "'"),
       cex.sub = 0.6)

  axis(1, at = 1:n_levels, labels = level_names, cex.axis = 0.8)
  axis(2, las = 1, at = round(seq(0, max(upper_bounds, na.rm = TRUE) * 1.2, by = max(upper_bounds, na.rm = TRUE) * 0.2), 2), cex.axis = 0.8)

  for (i in 1:length(odds_ratios)) {
    if (!is.na(lower_bounds[i]) && !is.na(upper_bounds[i])) {
      lines(c(i, i), c(lower_bounds[i], upper_bounds[i]), lwd = 1)
      text(i, lower_bounds[i], labels = round(lower_bounds[i], 2), pos = 2, cex = 0.45, offset = 0.3)
      text(i, upper_bounds[i], labels = round(upper_bounds[i], 2), pos = 2, cex = 0.45, offset = 0.3)
    }
    points(i, odds_ratios[i], pch = 16)
    text(i, odds_ratios[i], labels = round(odds_ratios[i], 2), pos = 4, cex = 0.55)
  }

  # Add reference line at y = 1 and tick mark
  abline(h = 1, lty = 2, col = "gray")
}
