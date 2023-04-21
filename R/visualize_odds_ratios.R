#' Visualize Odds Ratios for a 2xk Contingency Table (internal function)
#'
#' This function creates a plot of odds ratios with 95% confidence intervals
#' for a 2xk contingency table.
#'
#' @param ctable A 2xk contingency table as a matrix or data frame with row and column names. The table must have dimensions of at least 2x3.
#' @param reference_level The index of the reference level for odds ratio calculations (default: 1). The user must select the column level to serve as the reference level.
#' @param row_category The index of the row category to be used in odds ratio calculations (1 or 2). The user must select the row level to which the calculation of the odds ratios make reference (default: 1).
#' @param or.alpha The significance level used for the confidence intervals (default: 0.05).
#'
#' @return A plot of odds ratios with 95% confidence intervals for a 2xk contingency table.
#'
#' @keywords internal
#'
#' @importFrom graphics axis lines text
#' @importFrom stats chisq.test
#'
#'
visualize_odds_ratios <- function(ctable, reference_level = 1, row_category = 1, or.alpha = 0.05) {
  # Convert data frame to matrix if needed
  if (is.data.frame(ctable)) {
    ctable <- as.matrix(ctable)
  }

  if (!any(dim(ctable) == 2) || ncol(ctable) < 3) {
    stop("The input cross-tab must have dimensions of at least 2x3.")
  }

  if (row_category != 1 && row_category != 2) {
    stop("Invalid row_category. Must be either 1 or 2.")
  }

  # Extract row_category name
  rc_name <- rownames(ctable)[row_category]

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
    if (i != reference_level) {
      if (row_category == 1) {
        odds_ratio <- (ctable[1, i] * ctable[2, reference_level]) / (ctable[1, reference_level] * ctable[2, i])
      } else {
        odds_ratio <- (ctable[2, i] * ctable[1, reference_level]) / (ctable[2, reference_level] * ctable[1, i])
      }

      se <- sqrt(sum(1 / ctable[c(row_category, 3 - row_category), c(i, reference_level)]))
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

  # Perform Chi-square test
  chisq_test <- chisq.test(ctable)

  # Create the visualization
  plot(NULL, xlim = c(0.5, n_levels + 0.5), ylim = c(0, max(upper_bounds, na.rm = TRUE) * 1.2),
       xlab = "", ylab = "Odds Ratios", xaxt = 'n', yaxt = 'n',
       main = "Odds Ratios with 95% Confidence Intervals", cex.main = 0.9,
       sub = paste0("Chi-square: ", round(chisq_test$statistic, 2),
                    ", df: ", chisq_test$parameter,
                    ", p-value: ", format.pval(chisq_test$p.value, digits = 3),
                    "\nThe ORs are relative to the '", rc_name, "' row level",
                    "\nReference column level: ", level_names[reference_level]),
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

