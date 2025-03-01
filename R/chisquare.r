#' R function for Chi-square, (N-1) Chi-square, and G-Square test of independence, power calculation, measures of association, and standardised/moment-corrected
#' standardised/adjusted standardised residuals, visualisation of odds ratio in 2xk tables (where k >= 2)
#'
#' @description The function performs the chi-square test (both in its original format and in the N-1 version) and the G-square test of independence
#' on the input contingency table. It also calculates the retrospective power of the traditional chi-square test and various measures of categorical association for tables of any size,
#' returns standardised, moment-corrected standardised, adjusted standardised residuals (with indication of their significance), Quetelet Index,
#' IJ association factor, adjusted standardised counts, and the chi-square-maximising version of the input table. It also calculates relative and absolute
#' contributions to the chi-square statistic. The p value associated to the chi-square statistic
#' is also computed via both a permutation- and a Monte Carlo-based method. The 95 percent confidence interval around those p values is also provided.
#' Nicely-formatted output tables are rendered. Optionally, in 2xk tables (where k >= 2), a plot of the odds ratios can be rendered.\cr
#' Visit this \href{https://drive.google.com/file/d/1WxRUOpKUGcHW8OwMJLS_uy5QmLaCplCa/view?usp=sharing}{LINK} to access the package's vignette.\cr
#'
#' @details The function produces the following \strong{measures of categorical associations}:
#'  \itemize{
#'   \item Phi (with indication of the magnitude of the effect size; only for 2x2 tables)
#'   \item Phi max (used to compute Phi corrected; only for 2x2 tables)
#'   \item Phi corrected (with indication of the magnitude of the effect size; only for 2x2 tables)
#'   \item Phi signed (with indication of the magnitude of the effect size; only for 2x2 tables)
#'   \item Yule's Q (only for 2x2 tables, includes 95perc confidence interval, p-value, and indication of the magnitude of the effect)
#'   \item Yule's Y (only for 2x2 tables, includes 95perc confidence interval, p-value and indication of the magnitude of the effect)
#'   \item Odds ratio (for 2x2 tables, includes 95perc confidence interval, p value, and indication of the magnitude of the effect)
#'   \item Independent odds ratios (for tables larger than 2x2)
#'   \item Adjusted contingency coefficient C (with indication of the magnitude of the effect)
#'   \item Cramer's V (with 95perc confidence interval; includes indication of the magnitude of the effect)
#'   \item Cramer's V max (used to compute V corrected; for tables of any size)
#'   \item Cramer's V corrected (with indication of the magnitude of the effect)
#'   \item Cramer's V standardised (with indication of the magnitude of the effect)
#'   \item Bias-corrected Cramer's V (with indication of the magnitude of the effect)
#'   \item Cohen's w (with indication of the magnitude of the effect)
#'   \item W coefficient (includes 95perc confidence interval and magnitude of the effect)
#'   \item Goodman-Kruskal's lambda (both asymmetric and symmetric)
#'   \item Corrected version of lambda (both asymmetric and symmetric)
#'   \item Goodman-Kruskal's tau (asymmetric)
#' }
#'
#' \strong{Retrospective Power of the Traditional Chi-Square Test}\cr
#' The function calculates the retrospective power of the traditional chi-square test, which is the probability of correctly rejecting the null
#' hypothesis when it is false. The power is determined by the observed chi-square statistic, the sample size,
#' and the degrees of freedom, without explicitly calculating an effect size, following the method described by Oyeyemi et al. 2010.
#'
#' The degrees of freedom are calculated as (number of rows - 1) * (number of columns - 1). The alpha level is set by default at 0.05
#' and can be customized using the \code{power.alpha} parameter. The power is then estimated using the non-centrality parameter based
#' on the observed chi-square statistic.
#'
#' The calculation involves determining the critical chi-squared value based on the alpha level and degrees of freedom, and then
#' computing the probability that the chi-squared distribution with the given degrees of freedom exceeds this critical value.
#'
#' The resulting power value indicates how likely the test is to detect an effect if one exists. A power value close to 1 suggests a
#' high probability of detecting a true effect, while a lower value indicates a higher risk of a Type II error. Typically, a power
#' value of 0.8 or higher is considered robust in most research contexts.
#'
#'
#' \strong{Suggestion of a suitable chi-square testing method}\cr
#' The first rendered table includes a suggestion for the applicable chi-squared test method,
#' derived from an internal analysis of the input contingency table. The decision logic used is as follows:\cr
#'
#' For 2x2 tables:\cr
#' - if the grand total is equal to or larger than 5 times the number of cells,
#'   the traditional Chi-Square test is suggested. Permutation or Monte Carlo
#'   methods can also be considered.\cr
#'
#' - if the grand total is smaller than 5 times the number of cells, the minimum expected count is checked:\cr
#' (A) if it is equal to or larger than 1, the (N-1)/N adjusted Chi-Square test is
#'     suggested, with an option for Permutation or Monte Carlo methods.\cr
#' (B) if it is less than 1, the Permutation or Monte Carlo method is recommended.\cr
#'
#' For tables larger than 2x2:\cr
#' - the logic is similar to that for 2x2 tables, with the same criteria for
#'   suggesting the traditional Chi-Square test, the (N-1)/N adjusted test,
#'   or the Permutation or Monte Carlo methods.\cr
#'
#' The rationale of a threshold for the applicability of the traditional chi-square test corresponding to
#' 5 times the number of cells is based on the following.\cr
#'
#' Literature indicates that the traditional chi-squared test's validity is not as fragile as once thought,
#' especially when considering the average expected frequency across all cells in the cross-tab, rather than
#' the minimum expected value in any single cell. An average expected frequency of at least 5 across all
#' cells of the input table should be sufficient for maintaining the chi-square test's reliability at the
#' 0.05 significance level.
#'
#' As a consequence, a table's grand total equal to or larger than 5 times the number of cells should ensure the applicability
#' of the traditional chi-square test (at alpha 0.05).
#'
#' See: Roscoe-Byars 1971; Greenwood-Nikulin 1996; Zar 2014; Alberti 2024.
#'
#' For the rationale of the use of the (N-1)/N adjusted version of the chi-square test,
#' and for the permutation and Monte Carlo method, see below.
#'
#'
#' \strong{Chi-square statistics adjusted using the (N-1)/N adjustment}\cr
#' The adjustment is done by multiplying the chi-square statistics by (N-1)/N, where N is the table grand total (sample size). The p-value
#' of the corrected statistic is calculated the regular way (i.e., using the same degrees of freedom as in the traditional test).
#' The correction seems particularly relevant for tables where N is smaller than 20 and where the expected frequencies are equal
#' or larger than 1. The corrected chi-square test proves more conservative when the sample size is small.
#' As N increases, the term (N-1)/N approaches 1, making the adjusted chi-square value virtually equivalent to the unadjusted value.\cr
#'
#' See: Upton 1982; Rhoades-Overall 1982; Campbel 2007; Richardson 2011; Alberti 2024.
#'
#'
#' \strong{Permutation-based and Monte Carlo p-value for the chi-square statistic}\cr
#' The p-value of the observed chi-square statistic is also calculated on the basis of both a permutation-based and a
#' Monte Carlo approach. In the first case, the dataset is permuted \emph{B} times (1000 by default), whereas in the second method
#' \emph{B} establishes the number of random tables generated under the null hypothesis of independence (1000 by default).\cr
#'
#' As for the permutation method, the function does the following internally:\cr
#' (1) Converts the input dataset to long format and expands to individual observations; \cr
#' (2) Calculates the observed chi-squared statistic; \cr
#' (3) Randomly shuffles (B times) the labels of the levels of one variable, and recalculates chi-squared statistic for each shuffled dataset;\cr
#' (4) Computes the p-value based on the distribution of permuted statistics (see below).\cr
#'
#' For the rationale of the permutation-based approach, see for instance Agresti et al 2022.\cr
#'
#' For the rationale of the Monte Carlo approach, see for instance the description in Beh-Lombardo 2014: 62-64.\cr
#'
#' Both simulated p-values are calculated as follows: \cr
#'
#' \eqn{sum (chistat.simulated >= chisq.stat) / B}, where\cr
#'
#' \emph{chistat.simulated} is a vector storing the B chi-squared statistics generated under the Null Hypothesis, and\cr
#' \emph{chisq.stat} is the observed chi-squared statistic.\cr
#'
#' Both distributions can be optionally plotted setting the \code{graph} parameter to \code{TRUE}.
#'
#'
#' \strong{Confidence interval around the permutation-based and Monte Carlo p-value}\cr
#' The function calculates the 95 percent Confidence Interval around the simulated p-values.
#' The Wald CI quantifies the uncertainty around the simulated p-value estimate. For a 95 percent CI,
#' the standard z-value of 1.96 is used. The standard error for the estimated p-value is computed as the square root of
#' (estimated p-value * (1 - estimated p-value) / number of simulations-1).
#'
#' The lower and upper bounds of the CI are then calculated as follows:\cr
#' Lower Confidence Interval = estimated p-value - (z-value * standard error)\cr
#' Upper Confidence Interval = estimated p-value + (z-value * standard error)\cr
#'
#' Finally, the lower and upper CIs are clipped to lie within 0 and 1.
#'
#' The implemented procedure aligns with the one described at this link:
#' https://blogs.sas.com/content/iml/2015/10/28/simulation-exact-tables.html
#'
#' For the Wald confidence interval see Fagerland et al 2017.
#'
#'
#' \strong{Cells' relative contribution (in percent) to the chi-square statistic}\cr
#' The cells' relative contribution (in percent) to the chi-square statistic is calculated as:\cr
#'
#' \eqn{chisq.values / chisq.stat * 100}, where\cr
#'
#' \emph{chisq.values} and \emph{chisq.stat} are the chi-square
#' value in each individual cell of the table and the value of the chi-square statistic, respectively. The
#' \emph{average contribution} is calculated as \eqn{100 / (nr*nc)}, where \emph{nr} and \emph{nc} are the
#' number of rows and columns in the table respectively.\cr
#'
#'
#' \strong{Cells' absolute contribution (in percent) to the chi-square statistic}\cr
#' The cells' absolute contribution (in percent) to the chi-square statistic is calculated as:\cr
#'
#' \eqn{chisq.values / n * 100}, where\cr
#'
#' \emph{chisq.values} and \emph{n} are the chi-square
#' value in each individual cell of the table and the table's grant total, respectively. The
#' \emph{average contribution} is calculated as sum of all the absolute contributions divided by the number of cells in
#' the table.\cr
#'
#' For both the relative and absolute contributions to the chi-square, see: Beasley-Schumacker 1995: 90.
#'
#'
#' \strong{Chi-square-maximising table}\cr
#' The chi-square-maximising table is the version of the input cross-tab that, while preserving the marginal
#' configuration, produces the largest divergence between the observed and the expected counts and, therefore,
#' the maximum chi-squared value achievable given the observed marginals.
#'
#' The table is worked out using the routine described by Berry and colleagues. This allocation routine effectively
#' maximises the chi-square statistic by concentrating (to the maximum extent possible given the marginals) the levels of
#' one variable into specific levels of the other.
#'
#' As Berry and colleagues have noted, there can be alternative
#' positions for the zeros in the chi-square-maximising table, but the obtained non-zero counts are the only
#' ones that allow the maximisation of the chi-squared statistic.
#'
#' The chi-square-maximising table is used to compute Cramer's V max, which is in turn used to compute
#' Cramer's V corrected (equal to phi-corrected for 2x2 tables; see below).
#'
#' On the chi-square-maximising table, see Berry et al. 2018.
#'
#'
#' \strong{Moment-corrected standardised residuals}\cr
#' The moment-corrected standardised residuals are calculated as follows: \cr
#'
#' \eqn{stand.res / (sqrt((nr-1)*(nc-1)/(nr*nc)))}, where\cr
#'
#' \emph{stand.res} is each cell's standardised residual, \emph{nr} and
#' \emph{nc} are the number of rows and columns respectively.\cr
#'
#' See Garcia-Perez-Nunez-Anton 2003: 827.\cr
#'
#'
#' \strong{Adjusted standardised residuals}\cr
#' The adjusted standardised residuals are calculated as follows: \cr
#'
#' \eqn{stand.res[i,j] / sqrt((1-sr[i]/n)*(1-sc[j]/n))}, where\cr
#'
#' \emph{stand.res} is the standardised residual for cell \emph{ij},
#' \emph{sr} is the row sum for row \emph{i}, \emph{sc} is the column sum for column \emph{j}, and
#' \emph{n} is the table grand total. The \emph{adjusted standardised residuals} should be used in place of
#' the standardised residuals since the latter are not truly standardised because they have a nonunit variance. The
#' standardised residuals therefore underestimate the divergence between the observed and the expected counts. The adjusted
#' standardised residuals (and the moment-corrected ones) correct that deficiency.\cr
#'
#' For more info see: Haberman 1973.
#'
#'
#' \strong{Significance of the residuals}\cr
#' The significance of the residuals (standardised, moment-corrected standardised, and adjusted standardised) is assessed using alpha 0.05 or, optionally
#' (by setting the parameter \code{adj.alpha} to \code{TRUE}),
#' using an adjusted alpha calculated using the Sidak's method:\cr
#'
#' \eqn{alpha.adj = 1-(1 - 0.05)^(1/(nr*nc))}, where\cr
#'
#' \emph{nr} and \emph{nc} are the number of rows and columns in the table respectively. The adjusted
#' alpha is then converted into a critical two-tailed z value. \cr
#'
#' See: Beasley-Schumacker 1995: 86, 89.
#'
#'
#' \strong{Quetelet Index and IJ association factor}\cr
#' The Quetelet Index expresses the relative change in probability in one variable when considering
#' the association with the other. The sign indicates the direction of the change, that is, whether the probability
#' increases or decreases, compared to the probability expected under the hypothesis of independence. The IJ association factor
#' is centred around 1.0 and represents the factor by which the probability in one variable changes when considering the association with the other.
#' A decrease in probability is indicated by a factor smaller than 1.\cr
#'
#' The Quetelet index is computes as:\cr
#'
#' \eqn{(observed freq / expected freq) - 1}\cr
#'
#' The IJ association factor is computed as:\cr
#'
#' \eqn{observed freq / expected freq}\cr
#'
#' The thresholds for an IJ factor indicating a noteworthy change in probability, based on Agresti 2013, are: "larger than 2.0" and
#' "smaller than 0.5". The thresholds for the Quetelet Index are based on those, and translate into "larger than 1.0" and "smaller than -0.50".
#' For example, a Quetelet index greater than 1.0 indicates that the observed probability is more than double the expected probability under independence,
#' corresponding to an increase of over 100 percent. Similarly, a Quetelet index less than -0.5 indicates that the observed probability is less than half the
#' expected probability, corresponding to a substantial decrease of over 50 percent.
#'
#' For the Quetelet Index, see: Mirkin 2023.
#'
#' For the IJ association factor, see: Agresti 2013; Fagerland et al 2017. Note that
#' the IJ association factor is called 'Pearson ratio' by Goodman 1996.
#'
#'
#' \strong{Adjusted standardised counts}\cr
#' The function computes adjusted standardised counts for the input contingency table. It first standardises the counts via
#' Iterative Proportional Fitting (see below) so that all row and column totals equal unity; then adjusts these
#' standardised counts by subtracting the table's grand mean (that is, the grand total divided by the number of cells).
#'
#' The resulting adjusted standardised counts isolate the core association structure of the contingency table by removing
#' the confounding effects of marginal distributions. This allows for meaningful comparisons across different
#' tables, even when the original tables have different marginal totals.
#'
#' In the adjusted standardised counts, a value of 0 indicates independence between the row and column variables for that cell.
#' Positive values indicate a higher association than expected under independence, while negative values
#' indicate a lower association than expected.
#'
#' Unlike standardised residuals, these adjusted standardised counts can be directly compared across tables of
#' the same size, making them particularly useful for comparative analyses across different samples, time periods, or contexts.
#'
#' It goes without saying that, in the process of table standardisation, any target marginals can be chosen.
#' The choice of target values affects the scale of the standardised counts but not their relative relationships.
#' By setting all the row and columns sums to unity, the function follows Rosenthal-Rosnow 2008.
#'
#' See: Fienberg 1971; Rosenthal-Rosnow 2008.
#'
#'
#' \strong{Phi corrected}\cr
#' To refine the phi coefficient, scholars have introduced a corrected version. It accounts for the fact that the original coefficient (1)
#' does not always have a maximum achievable value of unity since it depends on the marginal configuration, and therefore (2) it is not directly
#' comparable across tables with different marginals. To calculate phi-corrected, one first computes phi-max, which represents the
#' maximum possible value of phi under the given marginal totals. Phi-corrected is equal to phi/phi-max.
#'
#' For more details see: Cureton 1959; Liu 1980; Davenport et al. 1991; Rash et al. 2011; Alberti 2024.
#'
#' See also \emph{Chi-square-maximising table} above.
#'
#'
#' \strong{95perc confidence interval around Cramer's V}\cr
#' The calculation of the 95perc confidence interval around Cramer's V is based on Smithson 2003: 39-41, and builds on the R code made
#' available by the author on the web (http://www.michaelsmithson.online/stats/CIstuff/CI.html).
#'
#'
#' \strong{Table standardisation via Iterative Proportional Fitting and Cramer's V standardised}\cr
#' This version of the coefficient is computed on the standardised table, which is returned and rendered by the function.
#' The \emph{standardised table} is internally obtained via
#' the \emph{Iterative Proportional Fitting} routine described, for instance, in Reynold 1977: 32-33. Since a number of association measures,
#' among which Cramer's V, are affected by the skewness in the marginal distributions, the original table is first standardised and
#' then Cramer's V is computed.
#'
#' The rationale of the use of standardised tables as basis to compute Cramer's V is that coefficients calculated on standardised tables
#' are comparable across tables because the impact of different marginal distributions is controlled for.
#'
#' The value obtained by subtracting the ratio Cramer's V to Cramer's V standardised from 1 gives an idea of the reduction of the magnitude of
#' V due to the skewness of the marginal sums (multiplied by 100 can be interpreted as percentage reduction).
#'
#' The standardisation is performed so that the rows feature the same sums, and the columns
#' features the same sum, while keeping the table's grand total unchanged. This removes the effect of skewed marginal distributions,
#' while preserving the association structure of the original table.
#'
#' The target row and column marginals used in the standardisation process are set using the \code{marginal.type}, \code{custom.row.totals}, and
#' \code{custom.col.totals} parameters.
#'
#' In the iterative procedure, the convergence to the established marginals is reached when the counts obtained
#' in a given iteration do not differ from the ones obtained in the previous iteration by more than a threshold internally set to 0.001.
#' The maximum number of iterations is internally set to 10,000 to prevent infinite loop; after that, the convergence is deemed as failed.
#'
#' For table standardisation as a preliminary step towards coefficients comparison across tables, see
#' for example Smith 1976; Reynolds 1977; Liu 1980.
#'
#'
#' \strong{Bias-corrected Cramer's V}\cr
#' The bias-corrected Cramer's V is based on Bergsma 2013: 323–328.\cr
#'
#'
#' \strong{W coefficient}\cr
#' It addresses some limitations of Cramer's V. When the marginal probabilities are unevenly distributed, V may overstate the
#' strength of the association, proving pretty high even when the overall association is weak. W is based on the distance between observed
#' and expected frequencies. It uses the squared distance to adjust for the unevenness of the marginal distributions in the table.
#' The indication of the magnitude of the association is based on Cohen 1988 (see above).
#' Unlike Kvalseth 2018a, the calculation of the 95 percent confidence interval is based on a bootstrap approach (employing 10k resampled tables, and the 2.5th and 97.5th
#' percentiles of the bootstrap distribution).
#'
#' For more details see: Kvalseth 2018a.
#'
#'
#' \strong{Indication of the magnitude of the association as indicated by the chi-squared-based coefficients}\cr
#' The function provides indication of the magnitude of the association (effect size) for the phi, phi signed, phi corrected, Cadj, Cramer's V,
#' V corrected, V standardised, V bias-corrected, Cohen's w, and W.
#'
#' The verbal articulation of the effect size is based on Cohen 1988.
#'
#' Phi, phi signed, phi corrected and w are assessed against the well-known Cohen's classification
#' scheme's thresholds (small 0.1, medium 0.3, large 0.5). For input cross-tabs larger than 2x2, the Cadj, Cramer's V, V corrected, V standardised, V bias-corrected,
#' and W coefficients are assessed against thresholds that depend on the table's df, which (as per Cohen 1988) correspond to the smaller between the
#' rows and columns number, minus 1. On the basis of the table's df, the three thresholds are calculated as follows:
#'
#' small effect: 0.100 / sqrt(min(nr,nc)-1)\cr
#' medium effect: 0.300 / sqrt(min(nr,nc)-1)\cr
#' large effect: 0.500 / sqrt(min(nr,nc)-1)\cr
#'
#' where nr and nc are the number of rows and number of columns respectively, and min(nr,nc)-1 corresponds to the table's df.
#' Essentially, the thresholds for a small, medium, and large effect are computed by dividing the Cohen's thresholds for a 2x2 table (df=1)
#' by the square root of the input table's df.
#'
#' Consider a V value of (say) 0.350; its effect size interpretation changes based on the table's dimension:\cr
#'
#' for a 2x2 table, 0.350 corresponds to a "medium" effect;\cr
#' for a 3x3 table, 0.350 still corresponds to a "medium" effect;\cr
#' for a 4x4 table, 0.350 corresponds to a "large" effect.
#'
#' The examples illustrate that for the same (say) V value, the interpreted effect size can shift from "medium" in a smaller table to "large" in a larger table.
#' In simpler terms, the threshold for determining a "large" effect, for instance, becomes more accessible to reach as the table's size increases.
#'
#' It is crucial to be aware of this as it highlights that the same coefficient value can imply different magnitudes of effect depending on the table's size.
#'
#' See: Cohen 1988; Sheskin 2011; Alberti 2024.
#'
#'
#' \strong{Corrected Goodman-Kruskal's lambda}\cr
#' The corrected Goodman-Kruskal's lambda adeptly addresses skewed or unbalanced marginal probabilities which create problems to the traditional lambda.
#' By emphasizing categories with higher probabilities through a process of squaring maximum probabilities and normalizing with marginal probabilities, this refined
#' coefficient addresses inherent limitations of lambda.
#'
#' For more details see: Kvalseth 2018b.
#'
#'
#' \strong{Magnitude of the association as indicated by Yule's Q and Yule's Y}\cr
#' Given the relationship between Q and Y and the odds ratio, the magnitude of the association indicated by Q and Y is based on the thresholds proposed by Ferguson 2009 for the
#' odds ratio (see below). Specifically, the thresholds for Q (in terms of absolute value) are:
#'
#' \itemize{
#'   \item |Q| < 0.330 - Negligible
#'   \item 0.330 <= |Q| < 0.500 - Small
#'   \item 0.500 <= |Q| < 0.600 - Medium
#'   \item |Q| >= 0.600 - Large
#' }
#'
#' For Yule's Y (absolute value):
#'
#'  \itemize{
#'   \item |Y| < 0.171 - Negligible
#'   \item 0.171 <= |Y| < 0.268 - Small
#'   \item 0.268 <= |Y| < 0.333 - Medium
#'   \item |Y| >= 0.333 - Large
#' }
#'
#' Yule's Q has a probabilistic interpretation. It tells us how much more likely is to draw pairs of individuals who share the same characteristics
#' (concordant pairs) as opposed to drawing pairs who do not (discordant pairs).
#'
#' Yule's Y represents the difference between the probabilities in the diagonal and off-diagonal cells, in the equivalent
#' symmetric tables with equal marginals. In other words, Y measures the extent to which the probability of falling in the diagonal cells exceeds the probability
#' of falling in the off-diagonal cells, in the standardised version of the table.
#'
#' Note that, if either cross-tab's diagonal features 0, the \emph{Haldane-Anscombe correction} is applied to both Q and Y.
#' This results in a coefficient approaching, but not reaching, 1 in such circumstances, with the exact value depending
#' on the table's frequencies. For more info on the correction, see the \emph{Odds Ratio} section below.
#'
#' On Yule's Q and Y, see Yule 1912.
#'
#' On Yule's Y, see also Reynolds 1977.
#'
#' On the probabilistic interpretation of Q, see: Davis 1971; Goodman-Kruskal 1979.
#'
#' On the sensitivity of Yule's Q and Pearson's phi to different types of associations, see Alberti 2024.
#'
#'
#' \strong{Odds Ratio}\cr
#' For 2x2 tables, a single odds ratio is computed.
#'
#' For tables larger than 2x2, independent odds ratios are computed for adjacent rows and columns, so representing
#' the minimal set of ORs needed to describe all pairwise relationships in the contingency table.
#'
#' For tables of size 2xk (where k >= 2), pairwise odds ratios can be plotted (along with their confidence interval) by
#' setting the \code{plor.or} parameter to \code{TRUE} (see the \emph{Odd Ratios plot} section further down).
#'
#' In all three scenarios, the \emph{Haldane-Anscombe correction} is applied when zeros are present along any of the table's diagonals.
#' This correction consists of adding 0.5 to every cell of the relevant 2x2 subtable before calculating the odds ratio. This ensures:
#'
#' - For 2x2 tables: The correction is applied to the entire table if either diagonal contains a zero.\cr
#' - For larger tables: The correction is applied to each 2x2 subtable used to calculate an independent odds ratio, if that subtable has a zero in either diagonal.\cr
#' - For 2xk tables: The correction is applied to each pairwise comparison that involves a zero in either diagonal.
#'
#' Note that, the \emph{Haldane-Anscombe correction} results in a finite (rather than infinite) odds ratio, with the exact value depending on the table's frequencies.
#'
#' On the Haldane-Anscombe correction, see: Fleiss et al 2003.
#'
#'
#' \strong{Odds Ratio effect size magnitude}\cr
#' The magnitude of the association indicated by the odds ratio is based on the thresholds (and corresponding reciprocals)
#' suggested by Ferguson 2009 (for other thresholds, see for instance Chen et al 2010):\cr
#'
#'  \itemize{
#'   \item OR < 2.0 - Negligible
#'   \item 2.0 <= OR < 3.0 - Small
#'   \item 3.0 <= OR < 4.0 - Medium
#'   \item OR >= 4.0 - Large
#' }
#'
#' As noted earlier, the effect size magnitude for Yule's Q and Yule's Y is based on the above odds ratio's thresholds.
#'
#'
#' \strong{Odd Ratios plot}\cr
#' For 2xk table, where k >= 2:\cr
#' by setting the \code{plor.or} parameter to \code{TRUE}, a plot showing the odds ratios and their 95percent confidence interval will be rendered.
#' The confidence level can be modified via the \code{or.alpha} parameter.
#'
#' The odds ratios are calculated for the column levels, and one of them
#' is to be selected by the user as a reference for comparison via the \code{reference.level} parameter (set to 1 by default).
#' Also, the user may want to select the row category to which the calculation of the odds ratios makes reference (using the \code{row.level} parameter,
#' which is set to 1 by default).
#'
#' As mentioned earlier, if any of the 2x2 subtables on which the odds ratio is calculated
#' features zeros along any of the diagonal, the \emph{Haldane-Anscombe} correction is applied.
#'
#' To better understand the rationale of plotting the odds ratios, consider the following example, which uses on the famous Titanic data:
#'
#' Create a 2x3 contingency table:\cr
#' \code{mytable <- matrix(c(123, 158, 528, 200, 119, 181), nrow = 2, byrow = TRUE)} \cr
#' \code{colnames(mytable) <- c("1st", "2nd", "3rd")} \cr
#' \code{rownames(mytable) <- c("Died", "Survived")} \cr
#'
#' Now, we perform the test and visualise the odds ratios:\cr
#' \code{chisquare(mytable, plot.or=TRUE, reference.level=1, row.level=1)} \cr
#'
#' In the rendered plot, we can see the odds ratios and confidence intervals for the second and third column level
#' (i.e., 2nd class and 3rd class) because the first column level has been selected as reference level. The odds ratios are calculated
#' making reference to the first row category (i.e., \emph{Died}). From the plot, we can see that, compared to the 1st class,
#' passengers on the 2nd class have 2.16 times larger odds of dying; passengers on the 3rd class have 4.74 times larger odds of dying
#' compared to the 1st class.
#'
#' Note that if we set the \code{row.level} parameter to \code{2}, we make reference to the second row category, i.e. \emph{Survived}:\cr
#' \code{chisquare(mytable, plot.or=TRUE, reference.level=1, row.level=2)}
#'
#' In the plot, we can see that passengers in the 2nd class have 0.46 times the odds of surviving of passengers in the 1st class, while
#' passengers from the 3rd class have 0.21 times the odds of surviving of those travelling in the 1st class.
#'
#'
#' \strong{Other measures of categorical association}\cr
#' For the other measures of categorical association provided by the function, see for example Sheskin 2011: 1415-1427.
#'
#'
#' \strong{Additional notes on calculations}:
#' \itemize{
#'   \item{the \strong{Phi} coefficient is based on the chi-square statistic as per Sheskin 2011's equation 16.21, whereas the
#'   \strong{Phi signed} is after Sheskin's equation 16.20;}
#'
#'   \item{the \strong{2-sided p value of Yule's Q} is calculated following Sheskin 2011's equation 16.24;}
#'
#'   \item{\strong{Cohen's w} is calculated as \eqn{V * sqrt(min(nr, nc)-1)}, where \emph{V} is Cramer's V, and \emph{nr} and \emph{nc}
#'   are the number of rows and columns respectively; see Sheskin 2011: 679;}
#'
#'   \item{the \strong{symmetric} version of \strong{Goodman-Kruskal's lambda} is calculated
#'   as per Reynolds 1984: 55-57;}
#'
#'   \item{\strong{Goodman-Kruskal's tau} is calculated as per Reynolds 1984: 57-60.}
#'}
#'
#'
#' @param data Dataframe containing the input contingency table.
#' @param B  Number of simulated tables to be used to calculate the permutation- and the Monte Carlo-based p value (1000 by default).
#' @param plot.or Takes TRUE or FALSE (default) if the user wants a plot of the odds ratios to be rendered (only for 2xk tables, where k >= 2).
#' @param reference.level The index of the column reference level for odds ratio calculations (default: 1).
#' The user must select the column level to serve as the reference level (only for 2xk tables, where k >= 2).
#' @param row.level The index of the row category to be used in odds ratio calculations (1 or 2; default: 1).
#' The user must select the row level to which the calculation of the odds ratios make reference (only for 2xk tables, where k >= 2).
#' @param or.alpha The significance level used for the odds ratios' confidence intervals (default: 0.05).
#' @param power.alpha The significance level used for the calculation of the power of the traditional chi-square test (default: 0.05).
#' @param adj.alpha  Takes TRUE or FALSE (default) if the user wants or does not want the significance level of the
#' residuals (standardised, adjusted standardised, and moment-corrected) to be corrected using the Sidak's adjustment method (see Details).
#' @param marginal.type Defines the target marginal sums used for table standardisation via Iterative Proportional Fitting.
#' It takes \emph{average} (default) to have target row and column marginals equal to the table's grand total divided
#' by the number of rows and columns, respectively; it takes \emph{percent} to have target marginals equal to fractions of a grand total set to 100.
#' @param custom.row.totals A vector of numbers indicating the target row marginals to be used for table standardisation via Iterative Proportional Fitting (NULL by default).
#' @param custom.col.totals A vector of numbers indicating the target column marginals to be used for table standardisation via Iterative Proportional Fitting(NULL by default).
#' @param format Takes \emph{short} (default) if the dataset is a dataframe storing a contingency table; if the
#' input dataset is a dataframe storing two columns that list the levels of the two categorical variables,
#' \emph{long} will preliminarily cross-tabulate the levels of the categorical variable in the 1st column against
#' the levels of the variable stored in the 2nd column.
#' @param graph Takes TRUE or FALSE (default) if the user wants or does not want to plot the permutation and Monte Carlo
#' distribution of the chi-square statistic across the number of simulated tables set by the B parameter.
#' @param oneplot Takes TRUE (default) or FALSE if the user wants or does not want to render of the permutation and Monte Carlo
#' distribution in the same plot.
#' @param tfs Numerical value to set the size of the font used in the main body of the various output tables (13 by default).
#'
#'
#' @return The function produces \strong{optional charts} (distribution of the permuted chi-square statistic
#' and a plot of the odds ratios between a reference column level and the other ones, the latter only for 2xk tables where k >= 2), and
#' a number of \strong{output tables} that are nicely formatted with the help of the \emph{gt} package.
#' The output tables are listed below:
#'
#'  \itemize{
#'   \item Input contingency table (with some essential analytical results annotated at the bottom)
#'   \item Expected frequencies
#'   \item Cells' chi-square value
#'   \item Cells' relative contribution (in percent) to the chi-square statistic (cells in RED feature a larger-than-average
#'   contribution)
#'   \item Cells' absolute contribution (in percent) to the chi-square statistic (colour same as above)
#'   \item Chi-square-maximising table (with indication of the associated chi-square value, that is, the maximum value
#'   of the chi-square statistic achievable given the table margins)
#'   \item Standardised residuals (RED for large significant residuals, BLUE for small significant residuals)
#'   \item Moment-corrected standardised residuals (colour same as above)
#'   \item Adjusted standardised residuals (colour same as above)
#'   \item Table of independent odds ratios (for tables larger than 2x2)
#'   \item Quetelet Indices
#'   \item IJ association factors
#'   \item Adjusted standardised counts
#'   \item Input contingency table standardised via Iterative Proportional Fitting
#'   \item Table of output statistics, p values, and association measures
#' }
#'
#' Also, the function returns a \strong{list containing the following elements}:
#' \itemize{
#'   \item \strong{input.table}:
#'     \itemize{
#'       \item \emph{crosstab}: input contingency table.
#'     }
#'     \item \strong{chi.sq.maxim.table}:
#'     \itemize{
#'       \item \emph{chi.sq.maximising.table}: chi-square-maximising table.
#'     }
#'   \item \strong{standardised.table}:
#'     \itemize{
#'       \item \emph{standard.table}: standardised table on which Cramer's V standardised is computed.
#'     }
#'   \item \strong{chi.sq.related.results}:
#'     \itemize{
#'       \item \emph{exp.freq}: table of expected frequencies.
#'       \item \emph{smallest.exp.freq}: smallest expected frequency.
#'       \item \emph{avrg.exp.freq}: average expected frequency.
#'       \item \emph{chisq.values}: cells' chi-square value.
#'       \item \emph{chisq.relat.contrib}: cells' relative contribution (in percent) to the chi-square statistic.
#'       \item \emph{chisq.abs.contrib}: cells' absolute contribution (in percent) to the chi-square statistic.
#'       \item \emph{chisq.statistic}: observed chi-square value.
#'       \item \emph{chisq.p.value}: p value of the chi-square statistic.
#'       \item \emph{chisq.max}: chi-square value computed on the chi-square-maximising table.
#'       \item \emph{chi.sq.power}: retrospective power of the traditional chi-square test.
#'       \item \emph{chisq.adj}: chi-square statistic adjusted using the (N-1)/N correction.
#'       \item \emph{chisq.adj.p.value}: p value of the adjusted chi-square statistic.
#'       \item \emph{chisq.p.value.perm}: permutation-based p value, based on B permuted tables.
#'       \item \emph{chisq.p.value.perm CI lower boundary}: lower boundary of the 95 percent CI around the permutation-based p value.
#'       \item \emph{chisq.p.value.perm CI upper boundary}: upper boundary of the 95 percent CI around the permutation-based p value.
#'       \item \emph{chisq.p.value.MC}: Monte Carlo p value, based on B random tables.
#'       \item \emph{chisq.p.value.MC CI lower boundary}: lower boundary of the 95 percent CI around the Monte Carlo p value.
#'       \item \emph{chisq.p.value.MC CI upper boundary}: upper boundary of the 95 percent CI around the Monte Carlo p value.
#'     }
#'   \item \strong{G.square}:
#'     \itemize{
#'       \item \emph{Gsq.statistic}: observed G-square value.
#'       \item \emph{Gsq.p.value}: p value of the G-square statistic.
#'     }
#'   \item \strong{post.hoc}:
#'     \itemize{
#'       \item \emph{stand.resid}: table of chi-square standardised residuals.
#'       \item \emph{mom.corr.stand.resid}: table of moment-corrected standardised residuals.
#'       \item \emph{adj.stand.resid}: table of adjusted standardised residuals.
#'       \item \emph{Quetelet.Index}: table of Quetelet indices.
#'       \item \emph{IJ.assoc.fact.}: table of IJ association factors.
#'       \item \emph{adj.stand.counts}: table of adjusted standardised counts.
#'     }
#'     \item \strong{margin.free.assoc.measures}:
#'     \itemize{
#'       \item \emph{Yule's Q}: Q coefficient (only for 2x2 tables).
#'       \item \emph{Yule's Q CI lower boundary}: lower boundary of the 95perc CI.
#'       \item \emph{Yule's Q CI upper boundary}: upper boundary of the 95perc CI.
#'       \item \emph{Yule's Q p.value}: 2-tailed p value of Yule's Q.
#'       \item \emph{Yule's Y}: Y coefficient (only for 2x2 tables).
#'       \item \emph{Yule's Y CI lower boundary}: lower boundary of the 95perc CI.
#'       \item \emph{Yule's Y CI upper boundary}: upper boundary of the 95perc CI.
#'       \item \emph{Yule's Y p.value}: 2-tailed p value of Yule's Y.
#'       \item \emph{Odds ratio}: odds ratio (for 2x2 tables).
#'       \item \emph{Odds ratio CI lower boundary}: lower boundary of the 95perc CI.
#'       \item \emph{Odds ratio CI upper boundary}: upper boundary of the 95perc CI.
#'       \item \emph{Odds ratio p.value}: p value of the odds ratio.
#'       \item \emph{ORs}: table of independent odds ratios (for tables larger than 2x2).
#'     }
#'   \item \strong{chi.sq.based.assoc.measures}:
#'     \itemize{
#'       \item \emph{Phi.signed}: signed Phi coefficient (only for 2x2 tables).
#'       \item \emph{Phi}: Phi coefficient (only for 2x2 tables).
#'       \item \emph{Phi.max}: Maximum value of Phi given the marginals (only for 2x2 tables).
#'       \item \emph{Phi.corr}: corrected Phi coefficient (equal to Phi/Phi max; only for 2x2 tables).
#'       \item \emph{Cadj}: adjusted contingency coefficient C.
#'       \item \emph{Cramer's V}: Cramer's V coefficient.
#'       \item \emph{Cramer's V CI lower boundary}: lower boundary of the 95perc CI.
#'       \item \emph{Cramer's V CI upper boundary}: upper boundary of the 95perc CI.
#'       \item \emph{Cramer's V max}: Maximum value of Cramer's V given the marginals.
#'       \item \emph{Cramer's V corr}: corrected V coefficient (equal to V/Vmax; for tables of any size).
#'       \item \emph{Cramer's V standard.}: Cramer's V computed on the standardised table.
#'       \item \emph{1-(Cramer's V/V standard.)}: value indicating the reduction of the magnitude of V due to the skewness of the marginal sums.
#'       \item \emph{Cramer's Vbc}: bias-corrected Cramer's V coefficient.
#'       \item \emph{w}: Cohen's w.
#'       \item \emph{W}: W coefficient.
#'       \item \emph{W CI lower boundary}: lower boundary of the 95perc CI.
#'       \item \emph{W CI upper boundary}: upper boundary of the 95perc CI.
#'     }
#'   \item \strong{PRE.assoc.measures}:
#'     \itemize{
#'       \item \emph{lambda (rows dep.)}: Goodman-Kruskal's lambda coefficient (considering the rows being the dependent variable).
#'       \item \emph{lambda (cols dep.)}: Goodman-Kruskal's lambda coefficient (considering the columns being the dependent variable).
#'       \item \emph{lambda (symmetric)}: Goodman-Kruskal's symmetric lambda coefficient.
#'       \item \emph{lambda corrected (rows dep.)}: corrected version of the lambda coefficient (considering the rows being the dependent variable).
#'       \item \emph{lambda corrected (cols dep.)}: corrected version of the lambda coefficient (considering the columns being the dependent variable).
#'       \item \emph{lambda corrected (symmetric)}: corrected version of the symmetric lambda coefficient.
#'       \item \emph{tau (rows dep.)}: Goodman-Kruskal's tau coefficient (considering the rows being the dependent variable).
#'       \item \emph{tau (cols dep.)}: Goodman-Kruskal's tau coefficient (considering the columns being the dependent variable).
#'     }
#' }
#'
#' \strong{Note} that the \emph{p-values} returned in the above list are expressed in scientific notation, whereas the ones reported in the
#' output table featuring the tests' result and measures of association are reported as broken down into classes (e.g., <0.05, or <0.01, etc),
#' with the exception of the Monte Carlo p-value and its CI.\cr
#'
#' The \strong{following examples}, which use in-built datasets, can be run to familiarise with the function:\cr
#'
#' -perform the test on the in-built 'social_class' dataset:\cr
#' \code{result <- chisquare(social_class)} \cr
#'
#' -perform the test on a 2x2 subset of the 'diseases' dataset:\cr
#' \code{mytable <- diseases[3:4,1:2]} \cr
#' \code{result <- chisquare(mytable)} \cr
#'
#' -perform the test on a 2x2 subset of the 'safety' dataset:\cr
#' \code{mytable <- safety[c(4,1),c(1,6)]} \cr
#' \code{result <- chisquare(mytable)} \cr
#'
#' -build a toy dataset in 'long' format (gender vs. opinion about death sentence):\cr
#' \code{mytable <- data.frame(GENDER=c(rep("F", 360), rep("M", 340)),
#' OPINION=c(rep("oppose", 235),
#'          rep("favour", 125),
#'          rep("oppose", 160),
#'          rep("favour", 180)))}
#'
#' -perform the test specifying that the input table is in 'long' format:\cr
#' \code{result <- chisquare(mytable, format="long")} \cr
#'
#'
#' @keywords chiperm
#'
#' @references Agresti, A. (2013). Categorical Data Analysis (3rd ed.). Wiley. ISBN 9780470463635.
#'
#' @references Agresti, A., Franklin, C., & Klingenberg, B. (2022). Statistics: The Art and Science of Learning from Data, (5th ed.). Pearson Education.
#'
#' @references Alberti, G. (2024). From Data to Insights: A Beginner's Guide to Cross-Tabulation Analysis. Routledge - CRC Press.
#'
#' @references Beh E.J., Lombardo R. 2014. Correspondence Analysis: Theory, Practice and New Strategies, Chichester, Wiley.
#'
#' @references Beasley TM and Schumacker RE. 1995. Multiple Regression Approach to Analyzing Contingency Tables: Post Hoc and Planned Comparison Procedures.
#' The Journal of Experimental Education, 64(1).
#'
#' @references Bergsma, W. 2013. A bias correction for Cramér's V and Tschuprow's T. Journal of the Korean Statistical Society. 42 (3).
#'
#' @references Berry, K. J., Johnston, J. E., & Mielke, P. W., Jr. (2018). The Measurement of Association: A Permutation Statistical Approach. Springer.
#'
#' @references Campbell, I. (2007). Chi-squared and Fisher–Irwin tests of two-by-two tables with small sample recommendations.
#' In Statistics in Medicine (Vol. 26, Issue 19, pp. 3661–3675).
#'
#' @references Chen, H., Cohen, P., and Chen, S. (2010). How Big is a Big Odds Ratio? Interpreting the Magnitudes of Odds Ratios in Epidemiological Studies.
#' In Communications in Statistics - Simulation and Computation (Vol. 39, Issue 4, pp. 860–864).
#'
#' @references Cohen, J. 1988. Statistical power analysis for the behavioral sciences (2nd ed). Hillsdale, N.J: L. Erlbaum Associates.
#'
#' @references Cureton, E. E. (1959). Note on phi/phimax. In Psychometrika (Vol. 24, Issue 1, pp. 89–91).
#'
#' @references Davenport, E. C., Jr., & El-Sanhurry, N. A. (1991). Phi/Phimax: Review and Synthesis. In Educational and Psychological
#' Measurement (Vol. 51, Issue 4, pp. 821–828).
#'
#' @references Davis, J. A. (1971). Elementary Survey Analysis. Prentice Hall. ISBN 9780132605472.
#'
#' @references Fagerland, M. W., Lydersen, S., & Laake, P. (2017). Statistical Analysis of Contingency Tables. CRC Press. ISBN 9781466588172.
#'
#' @references Ferguson, C. J. (2009). An effect size primer: A guide for clinicians and researchers. Professional Psychology:
#' Research and Practice, 40(5), 532–538.
#'
#' @references Fienberg, S. E. (1971). A statistical technique for historians: Standardizing tables of counts. The Journal
#' of Interdisciplinary History, 1(2), 305-315.
#'
#' @references Fleiss, J. L., Levin, B., & Paik, M. C. 2003. Statistical Methods for Rates and Proportions (3rd ed.). Wiley.
#'
#' @references Garcia-Perez, MA, and Nunez-Anton, V. 2003. Cellwise Residual Analysis in Two-Way Contingency Tables. Educational and Psychological Measurement, 63(5).
#'
#' @references Goodman, L. A. (1996). A Single General Method for the Analysis of Cross-Classified Data: Reconciliation and Synthesis of Some Methods of Pearson,
#' Yule, and Fisher, and also Some Methods of Correspondence Analysis and Association Analysis. Journal of the American Statistical Association, 91(433), 408-428.
#'
#' @references Goodman, L. A., & Kruskal, W. H. (1979). Measures of Association for Cross Classifications. Springer-Verlag. ISBN 9780387904436.
#'
#' @references Greenwood, P. E., & Nikulin, M. S. (1996). A guide to chi-squared testing. John Wiley & Sons.
#'
#' @references Haberman, S. J. (1973). The Analysis of Residuals in Cross-Classified Tables. In Biometrics (Vol. 29, Issue 1, p. 205).
#'
#' @references Kvålseth, T. O. (2018a). An alternative to Cramér’s coefficient of association. In Communications in Statistics - Theory and Methods (Vol. 47, Issue 23, pp. 5662–5674).
#'
#' @references Kvålseth, T. O. (2018b). Measuring association between nominal categorical variables: an alternative to the Goodman–Kruskal lambda. In Journal of Applied Statistics
#' (Vol. 45, Issue 6, pp. 1118–1132).
#'
#' @references Liu, R (1980). A Note on Phi-Coefficient Comparison. In Research in Higher Education (Vol. 13, No. 1, pp. 3-8).
#'
#' @references Mirkin, B. (2023). A straightforward approach to chi-squared analysis of associations in contingency tables.
#' In E. J. Beh, R. Lombardo, & J. G. Clavel (Eds.), Analysis of Categorical Data from Historical
#' Perspectives (Behaviormetrics: Quantitative Approaches to Human Behavior, vol. 17). Springer.
#'
#' @references Oyeyemi, G. M., Adewara, A. A., Adebola, F. B., & Salau, S. I. (2010). On the Estimation of Power and Sample Size in Test of Independence.
#'  In Asian Journal of Mathematics and Statistics (Vol. 3, Issue 3, pp. 139–146).
#'
#' @references Rasch, D., Kubinger, K. D., & Yanagida, T. (2011). Statistics in Psychology Using R and SPSS. Wiley.
#'
#' @references Reynolds, H. T. 1977. The Analysis of Cross-Classifications. New York: Free Press.
#'
#' @references Reynolds, H. T. 1984. Analysis of Nominal Data (Quantitative Applications in the Social Sciences) (1st ed.). SAGE Publications.
#'
#' @references Rhoades, H. M., & Overall, J. E. (1982). A sample size correction for Pearson chi-square in 2×2 contingency tables. In Psychological Bulletin (Vol. 91, Issue 2, pp. 418–423).
#'
#' @references Richardson, J. T. E. (2011). The analysis of 2 × 2 contingency tables-Yet again. In Statistics in Medicine (Vol. 30, Issue 8, pp. 890–890).
#'
#' @references Rosenthal, R., & Rosnow, R. L. (2008). Essentials of Behavioral Research: Methods and Data Analysis (3rd ed.). McGraw-Hill Higher Education.
#'
#' @references Roscoe, J. T., & Byars, J. A. (1971). An Investigation of the Restraints with Respect to Sample Size Commonly Imposed on the Use of the Chi-Square Statistic.
#' Journal of the American Statistical Association, 66(336), 755–759.
#'
#' @references Sheskin, D. J. 2011. Handbook of Parametric and Nonparametric Statistical Procedures, Fifth Edition (5th ed.). Chapman and Hall/CRC.
#'
#' @references Smith, K. W. (1976). Marginal Standardization and Table Shrinking: Aids in the Traditional Analysis of Contingency Tables. Social Forces, 54(3), 669-693.
#'
#' @references Smithson M.J. 2003. Confidence Intervals, Quantitative Applications in the Social Sciences Series, No. 140. Thousand Oaks, CA: Sage.
#'
#' @references Upton, G. J. G. (1982). A Comparison of Alternative Tests for the 2 × 2 Comparative Trial. In Journal of the Royal Statistical Society.
#' Series A (General) (Vol. 145, Issue 1, p. 86).
#'
#' @references Yule, G. U. (1912). On the methods of measuring association between two attributes. Journal of the Royal Statistical Society, 75(6), 579–652.
#'
#' @references Zar, J. H. (2014). Biostatistical analysis (5th ed.). Pearson New International Edition.
#'
#' @export
#'
#' @importFrom gt gt cols_align tab_header md tab_source_note tab_options pct
#' @importFrom stats pchisq qchisq addmargins r2dtable pnorm quantile qnorm rmultinom
#' @importFrom graphics abline points hist rug layout par
#'
#'
#' @examples
#'
#' # Perform the analysis on a 2x2 subset of the in-built 'social_class' dataset
#' result <- chisquare(social_class[c(1:2), c(1:2)], B=9)
#'
#'
chisquare <- function(data, B = 1000, plot.or = FALSE, reference.level = 1, row.level = 1, or.alpha = 0.05,
                      power.alpha = 0.05, adj.alpha = FALSE, marginal.type = "average",
                      custom.row.totals = NULL, custom.col.totals = NULL, format = "short",
                      graph = FALSE, oneplot = TRUE, tfs = 13){
  VALUE <- NULL
  df <- data

  #if 'format' is equal to long
  if (format=="long") {
    #cross-tabulate the first column vs the second column
    from.long.to.short <- table(df[,1], df[,2])
    #turn the 'table' format to 'dataframe'
    df <- as.data.frame.matrix(from.long.to.short)
  }

  #n of rows
  nr <- as.integer(nrow(df))
  #n of columns
  nc <- as.integer(ncol(df))
  #n rows total
  sr <- rowSums(df)
  #n columns total
  sc <- colSums(df)
  #grand  total
  n <- sum(df)
  #degrees of freedom
  degrees.of.f <- (nr-1)*(nc-1)

  #average expected count
  avrg.expt.count <- round(n/(nr*nc),3)

  #if alpha.adj is TRUE, adjust the value of alpha using the Sidak's method,
  #and calculate the two-tailed critical value to be used as threshold for the
  #significance of the standardised and adjusted standardised residuals
  if (adj.alpha==TRUE) {
    alpha.adj = 1-(1 - 0.05)^(1/(nr*nc))
    z.crit.two.tailed <- qnorm(alpha.adj/2, mean = 0, sd = 1, lower.tail = F)

    #define a note to be used later on as annotation for the tables of residuals
    note.for.residuals <-  paste0("*BLUE: significant negative residuals (< ", round(0-z.crit.two.tailed,3), ") <br>
                                          RED: significant positive residuals (> ", round(z.crit.two.tailed,3), ") <br> <br>
                                          Highlighted residuals are significant at least at alpha ", round(alpha.adj,3), " (Sidak's alpha adjustment applied)*")
  } else{
    alpha.adj <- 0.05
    z.crit.two.tailed <- 1.96

    #define a note to be used later on as annotation for the tables of residuals
    note.for.residuals <-  paste0("*BLUE: significant negative residuals (< -1.96)<br>
                                          RED: significant positive residuals (> 1.96) <br> <br>
                                          Highlighted residuals are significant at least at alpha 0.05*")
  }

  #create a function to calculate the chi-sq statistic, to be also
  #used later on to calculate the simulated chi-sq statistic
  calc <- function(x){
    # Re-calculate the grand total based on the current table;
    # this is important because the function is also used to compute the chi-sq test
    # on the standardised table, and its grand total may not correspond to the
    # grand total of the original inputted table
    n <- sum(x)

    #expected frequencies
    exp.freq <- round(outer(rowSums(x), colSums(x), "*")/n, 3)

    #table with chi-square values per each cell
    chisq.values <- (x - exp.freq)^2/exp.freq

    #chi-square statistic
    chisq.stat <- round(sum(chisq.values),3)

    results <- list("chisq.stat"=chisq.stat,
                    "chisq.values"=chisq.values,
                    "exp.freq"=exp.freq)
    return(results)
  }

  #put the above function to work and extract the chi-sq statistic
  #the chi-sq values, and the expected frequencies, and the adjusted chi-sq to be used later on
  chisq.stat <- calc(df)$chisq.stat
  chisq.values <- round(calc(df)$chisq.values,3)
  exp.freq <- round(calc(df)$exp.freq,3)
  chisq.adj <- round(chisq.stat * ((n-1)/n), 3)

  #p values for the chi.sq and chi-sq adjusted statistics
  p <- as.numeric(format(pchisq(q=chisq.stat, df=degrees.of.f, lower.tail=FALSE), scientific = T))
  p.chisq.adj <- as.numeric(format(pchisq(q=chisq.adj, df=degrees.of.f, lower.tail=FALSE), scientific = T))

  p.to.report <- ifelse(p < 0.001, "< 0.001",
                        ifelse(p < 0.01, "< 0.01",
                               ifelse(p < 0.05, "< 0.05",
                                      round(p, 3))))

  p.chisq.adj.to.report <- ifelse(p.chisq.adj < 0.001, "< 0.001",
                                  ifelse(p.chisq.adj < 0.01, "< 0.01",
                                         ifelse(p.chisq.adj < 0.05, "< 0.05",
                                                round(p.chisq.adj, 3))))


  ##calculate the chi-sq statistic B times, using the above-defined 'calc' function##
  extract <- function(x)  calc(x)$chisq.stat
  chistat.mc <- sapply(r2dtable(B, sr, sc), extract)

  #calculate the p value of the observed chi-sq on the basis of the B Monte Carlo chi-sq statistics
  p.uppertail <- sum (chistat.mc >= chisq.stat) / B

  #p.uppertail.to.report <- ifelse(p.uppertail < 0.001, "< 0.001",
  #ifelse(p.uppertail < 0.01, "< 0.01",
  #ifelse(p.uppertail < 0.05, "< 0.05",
  #round(p.uppertail, 3))))

  # Calculation of the confidence interval around the Monte Carlo p-value
  SE <- sqrt(p.uppertail * (1 - p.uppertail) / (B - 1))

  # set the apha value for a 95% CI
  wald.aplha <- 0.05

  # Z-value for the given alpha level
  z_value <- qnorm(1 - wald.aplha / 2)

  # Confidence interval calculation
  lower_ci <- p.uppertail - z_value * SE
  upper_ci <- p.uppertail + z_value * SE

  # Ensure CI is within [0, 1]
  lower_ci <- max(lower_ci, 0)
  upper_ci <- min(upper_ci, 1)


  ##carry out the permutation-based chi-sq test##

  # Step 1: Convert input table to long format
  long.format <- as.data.frame(as.table(as.matrix(df)))
  names(long.format) <- c("Variable1", "Variable2", "Count")

  # Step 2: Expand to individual observations
  expanded.data <- long.format[rep(seq_len(nrow(long.format)), long.format$Count), 1:2]

  # Calculate the original chi-squared statistic, using the above-defined 'calc()' function
  original.chi.squared <- calc(df)$chisq.stat

  # Store permuted chi-squared values
  permuted.chi.squared <- numeric(B)

  for (i in 1:B) {
    # Shuffle the labels of one of the variables
    shuffled.data <- expanded.data
    shuffled.data$Variable2 <- sample(shuffled.data$Variable2)

    # Re-cross-tabulate and compute chi-squared
    shuffled.table <- table(shuffled.data)
    permuted.chi.squared[i] <- calc(shuffled.table)$chisq.stat
  }

  # Calculate the p-value
  p.value.permut <- sum(permuted.chi.squared >= original.chi.squared) / B

  # Calculation of the confidence interval around the permuted p-value
  SE.perm <- sqrt(p.value.permut * (1 - p.value.permut) / (B - 1))

  # set the apha value for a 95% CI
  wald.aplha.perm <- 0.05

  # Z-value for the given alpha level
  z_value.perm <- qnorm(1 - wald.aplha / 2)

  # Confidence interval calculation
  lower_ci.perm <- p.value.permut - z_value.perm * SE.perm
  upper_ci.perm <- p.value.permut + z_value.perm * SE.perm

  # Ensure CI is within [0, 1]
  lower_ci.perm <- max(lower_ci.perm, 0)
  upper_ci <- min(upper_ci.perm, 1)


  ##G-squared##
  #if any of the cells in the input table is equal to 0,
  #replace the 0s with a small non-zero value because otherwise
  #the log of 0 would be undefined
  if(any(df==0)){
    df[df == 0] <- 0.001
  }
  Gsquared <- round(2*sum(df*log(df/exp.freq)),3)

  p.Gsquared <- as.numeric(format(pchisq(q=Gsquared, df=degrees.of.f, lower.tail=FALSE), scientific = T))

  p.Gsquared.to.report <- ifelse(p.Gsquared < 0.001, "< 0.001",
                                 ifelse(p.Gsquared < 0.01, "< 0.01",
                                        ifelse(p.Gsquared < 0.05, "< 0.05",
                                               round(p.Gsquared,3))))

  ## To play safe, assign the input dataset to the object 'df' again,
  # because of the preceding addition of 0.001 in case of 0s
  df <- data


  ## Define a function to calculate different versions of Phi ##
  calculate_phi <- function(df) {
    # extract cell counts from the table
    a <- as.numeric(df[1,1])
    b <- as.numeric(df[1,2])
    c <- as.numeric(df[2,1])
    d <- as.numeric(df[2,2])

    # calculate the phi coefficient
    phi <- round(sqrt(chisq.stat / n), 3)

    # calculate the phi.signed
    phi_signed <- round(((a*d)-(b*c)) / sqrt((a+b)*(c+d)*(a+c)*(b+d)),3)

    # calculate phi_max using the internal compute_phi_max() function
    phi_max <- compute_phi_max(df)

    #calculate phi.corrected
    phi_corr <- round(phi / phi_max, 3)

    return(list(phi = phi, phi_corr = phi_corr, phi_max = phi_max, phi_signed=phi_signed))
  }

  # Function to determine the effect size based on the phi coefficient value
  get_effect_size_label <- function(phi_coefficient) {
    # Define breakpoints for effect sizes
    breakpoints_phi <- c(-Inf, 0.10, 0.30, 0.50, Inf)
    labels_phi <- c("(negligible effect)", "(small effect)", "(medium effect)", "(large effect)")

    # Use the cut() function to determine the effect size for phi_coefficient
    effect_size <- cut(phi_coefficient, breaks = breakpoints_phi, labels = labels_phi, right = TRUE)

    # Return the effect size
    return(as.character(effect_size))
  }

  # put the 'calculate_phi' function to work and get the different phis (only for 2x2 tables)
  if (nr == 2 & nc == 2) {
    phis <- calculate_phi(df)
    phi <- phis$phi
    phi.signed <- phis$phi_signed
    phi.max <-  round(phis$phi_max,3)
    phi.corr <- phis$phi_corr

    # Get effect size labels
    phi_effect_size <- get_effect_size_label(abs(phi))
    phi_signed_effect_size <- get_effect_size_label(abs(phi.signed))
    phi_corr_effect_size <- get_effect_size_label(abs(phi.corr))

  } else {
    phi <- "-"
    phi.signed <- "-"
    phi.max <- "-"
    phi.corr <- "-"
    # If phi coefficients aren't calculated, set their effect sizes to "-"
    phi_effect_size <- ""
    phi_signed_effect_size <- ""
    phi_corr_effect_size <- ""
  }


  ##Yule's Q##
  if (nr==2 & nc==2) {

    # Check for zeros along the diagonal of the 2x2 table
    if (df[1,1] * df[2,2] == 0 || df[1,2] * df[2,1] == 0) {
      # In case of zeros along any of the diagonals, add 0.5 to every cell
      df.to.use <- df + 0.5
    } else {
      df.to.use <- df
    }

    Q <- round(((df.to.use[1,1]*df.to.use[2,2]) - (df.to.use[1,2]*df.to.use[2,1])) /
                 ((df.to.use[1,1]*df.to.use[2,2]) + (df.to.use[1,2]*df.to.use[2,1])),3)

    Q.se <- sqrt((1/4)*(1-Q^2)^2*(1/df.to.use[1,1]+1/df.to.use[1,2]+1/df.to.use[2,1]+1/df.to.use[2,2]))

    Q.z <- Q / Q.se

    Q.p <- as.numeric(format(2*pnorm(q=abs(Q.z), lower.tail=FALSE), scientific=T))

    Q.p.to.report <- ifelse(Q.p < 0.001, "< 0.001",
                            ifelse(Q.p < 0.01, "< 0.01",
                                   ifelse(Q.p < 0.05, "< 0.05",
                                          round(Q.p, 3))))

    #calculate the lower and upper limit for the 95% CI around Q
    Q.lower <- round(Q - 1.96 * Q.se,3)
    Q.upper <- round(Q + 1.96 * Q.se,3)

    #create a vector to store the results to be later reported in the output table
    report.of.Q <- paste0(Q, " (95% CI ", Q.lower,"-", Q.upper, "; p: ", Q.p.to.report, ")")

  } else {

    Q <- "-"
    Q.lower <- "-"
    Q.upper <- "-"
    Q.p <- "-"
    report.of.Q <- "-"
  }

  if (nr == 2 & nc == 2) {
    # Define thresholds for Yule's Q
    Q_small <- 0.330
    Q_medium <- 0.500
    Q_large <- 0.600

    # Use absolute value of Q to focus on magnitude only
    abs_Q <- abs(Q)

    # Define breakpoints and labels for effect sizes
    breakpoints_Q <- c(-Inf, Q_small, Q_medium, Q_large, Inf)
    labels_Q <- c("(negligible effect)", "(small effect)", "(medium effect)", "(large effect)")

    # Use cut() function to categorize the effect size
    Q.effect.size <- cut(abs_Q, breaks = breakpoints_Q, labels = labels_Q, include.lowest = TRUE)
  } else {
    Q.effect.size <- ""
  }


  ##Yule's Y##
  if (nr==2 & nc==2) {

    # Check for zeros along the diagonal of the 2x2 table
    if (df[1,1] * df[2,2] == 0 || df[1,2] * df[2,1] == 0) {
      # In case of zeros along any of the diagonals, add 0.5 to every cell
      df.to.use <- df + 0.5
    } else {
      df.to.use <- df
    }

    Y <- round((sqrt(df.to.use[1,1]*df.to.use[2,2]) - sqrt(df.to.use[1,2]*df.to.use[2,1])) /
                 (sqrt(df.to.use[1,1]*df.to.use[2,2]) + sqrt(df.to.use[1,2]*df.to.use[2,1])), 3)

    Y.se <- sqrt((1/16)*(1-Y^2)^2*(1/df.to.use[1,1]+1/df.to.use[1,2]+1/df.to.use[2,1]+1/df.to.use[2,2]))

    Y.z <- Y / Y.se

    Y.p <- as.numeric(format(2*pnorm(q=abs(Y.z), lower.tail=FALSE), scientific=T))

    Y.p.to.report <- ifelse(Y.p < 0.001, "< 0.001",
                            ifelse(Y.p < 0.01, "< 0.01",
                                   ifelse(Y.p < 0.05, "< 0.05",
                                          round(Y.p, 3))))

    #calculate the lower and upper limit for the 95% CI around Y
    Y.lower <- round(Y - 1.96 * Y.se,3)
    Y.upper <- round(Y + 1.96 * Y.se,3)

    #create a vector to store the results to be later reported in the output table
    report.of.Y <- paste0(Y, " (95% CI ", Y.lower,"-", Y.upper, "; p: ", Y.p.to.report, ")")

  } else {

    Y <- "-"
    Y.lower <- "-"
    Y.upper <- "-"
    Y.p <- "-"
    report.of.Y <- "-"
  }

  if (nr == 2 & nc == 2) {
    # Define thresholds for Yule's Y
    Y_small <- 0.171
    Y_medium <- 0.268
    Y_large <- 0.333

    # Use absolute value of Y to focus on magnitude only
    abs_Y <- abs(Y)

    # Define breakpoints and labels for effect sizes
    breakpoints_Y <- c(-Inf, Y_small, Y_medium, Y_large, Inf)
    labels_Y <- c("(negligible effect)", "(small effect)", "(medium effect)", "(large effect)")

    # Use cut() function to categorize the effect size
    Y.effect.size <- cut(abs_Y, breaks = breakpoints_Y, labels = labels_Y, include.lowest = TRUE)
  } else {
    Y.effect.size <- ""
  }


  ##Define Cohen's thresholds to be used for C, V, Vbc, and W, on the basis of the table's df##
  Cohen_small <- 0.100 / sqrt(min(nr,nc)-1)
  Cohen_medium <- 0.300 / sqrt(min(nr,nc)-1)
  Cohen_large <- 0.500 / sqrt(min(nr,nc)-1)

  # Define breakpoints for effect sizes
  breakpoints_Cohen <- c(-Inf, Cohen_small, Cohen_medium, Cohen_large, Inf)
  labels_Cohen <- c("(negligible effect)", "(small effect)", "(medium effect)", "(large effect)")


  ##adjusted contingency coefficient##
  C <- sqrt(chisq.stat / (n + chisq.stat))
  Cmax <- sqrt((min(nr,nc)-1) / min(nr,nc))
  Cadj <- round(C / Cmax,3)

  # Use the cut() function to determine the effect size
  Cadj.effect.size <- cut(Cadj, breaks = breakpoints_Cohen, labels = labels_Cohen, right = FALSE)


  ##Cramer's V##
  V <- round(sqrt(chisq.stat / (n * min(nr-1, nc-1))), 3)

  # Use the cut() function to determine the effect size
  V.effect.size <- cut(V, breaks = breakpoints_Cohen, labels = labels_Cohen, right = FALSE)

  #95percent CI around Cramer's V
  #fuction to calculate Delta Lower (after Smithson)
  lochi <- function(chival,df,conf){
    ulim <- 1 - (1-conf)/2
    lc <- c(.001,chival/2,chival)
    while(pchisq(chival,df,lc[1])<ulim) {
      if(pchisq(chival,df)<ulim)
        return(c(0,pchisq(chival,df)))
      lc <- c(lc[1]/4,lc[1],lc[3])
    }
    diff <- 1
    while(diff > .00001) {
      if(pchisq(chival,df,lc[2])<ulim)
        lc <- c(lc[1],(lc[1]+lc[2])/2,lc[2])
      else lc <- c(lc[2],(lc[2]+lc[3])/2,lc[3])
      diff <- abs(pchisq(chival,df,lc[2]) - ulim)
      ucdf <- pchisq(chival,df,lc[2])
    }
    c(lc[2],ucdf)
  }
  #fuction to calculate Delta Upper (after Smithson)
  hichi <- function(chival,df,conf){
    uc <- c(chival,2*chival,3*chival)
    llim <- (1-conf)/2
    while(pchisq(chival,df,uc[1])<llim) {
      uc <- c(uc[1]/4,uc[1],uc[3])
    }
    while(pchisq(chival,df,uc[3])>llim) {
      uc <- c(uc[1],uc[3],uc[3]+chival)
    }
    diff <- 1
    while(diff > .00001) {
      if(pchisq(chival,df,uc[2])<llim)
        uc <- c(uc[1],(uc[1]+uc[2])/2,uc[2])
      else uc <- c(uc[2],(uc[2]+uc[3])/2,uc[3])
      diff <- abs(pchisq(chival,df,uc[2]) - llim)
      lcdf <- pchisq(chival,df,uc[2])
    }
    c(uc[2],lcdf)
  }

  #put the above functions to work to get delta.lower and delta.upper
  delta.lower <- lochi(chival=chisq.stat, df=degrees.of.f, conf=0.95)[1]
  delta.upper <- hichi(chival=chisq.stat, df=degrees.of.f, conf=0.95)[1]

  #calculate the lower and upper limit for the 95% CI around V (after Smithson)
  V.lower <- round(sqrt((delta.lower + degrees.of.f) / (n * (min(nr, nc)-1))),3)
  V.upper <- round(sqrt((delta.upper + degrees.of.f) / (n * (min(nr, nc)-1))),3)

  #create a vector to store the results to be later reported in the output table
  report.of.V <- paste0(V, " (95% CI ", V.lower,"-", V.upper, ")")


  ##compute Cramer's V on the standardised input table##
  # use the ancillary function to compute the standardised table
  standardisation.res <- standardize_table(df, marginal.type = marginal.type, custom.row.totals = custom.row.totals, custom.col.totals = custom.col.totals)
  table.stand <- standardisation.res$table.stand

  # extract the iterations number to be used later on in the 'gt' table annotation
  n.iterations <- standardisation.res$n.iterations

  # compute the chi.sq statistic on the stand table in order to compute V on that table
  chisq.stat.stand <- calc(table.stand)$chisq.stat

  # calculate V on the standardised table
  V.stand <- round(sqrt(chisq.stat.stand / (sum(table.stand) * min(nr-1, nc-1))), 3)

  # Use the cut() function to determine the effect size
  V.stand.effect.size <- cut(V.stand, breaks = breakpoints_Cohen, labels = labels_Cohen, right = FALSE)

  # compute the ratio between V and Vstand
  VtoV.stand <- round(1-(V/V.stand),3)


  ##bias-corrected Cramer's V##
  phi.new <- max(0, (chisq.stat / n) - ((nc-1)*(nr-1)/(n-1)))
  k.new  <-  nc-((nc-1)^2 / (n-1))
  r.new  <-  nr-((nr-1)^2 / (n-1))
  V.bc <- round(sqrt(phi.new / min(k.new-1, r.new-1)), 3)

  # Use the cut() function to determine the effect size
  Vbc.effect.size <- cut(V.bc, breaks = breakpoints_Cohen, labels = labels_Cohen, right = FALSE)


  ## Define a function to calculate the W coefficient ##
  calculate_W <- function(df) {

    # Extract joint probabilities
    pij <- df / sum(df)

    # Extract marginal probabilities
    pi_plus <- rowSums(pij)
    pj_plus <- colSums(pij)

    # Calculate d^2
    d2 <- sum((pij - outer(pi_plus, pj_plus))^2)

    # Calculate the normalization term
    normalization <- d2 - sum(pij^2) + min(sum(pi_plus^2), sum(pj_plus^2))

    # Calculate W
    W_hat <- sqrt(d2) / sqrt(normalization)
    return(W_hat)
  }

  # Put the above function to work and get the W coefficient
  W <- round(calculate_W(df),3)

  # Bootstrap the W coefficient
  Bboot <- 10000
  W_values <- numeric(Bboot)
  N <- sum(df)
  for (i in 1:Bboot) {
    resampled_vector <- rmultinom(1, N, prob = as.vector(as.matrix(df)/N))
    resampled_table <- matrix(resampled_vector, nrow = nrow(df))
    W_values[i] <- calculate_W(resampled_table)
  }

  # Calculate 95% CI
  CI <- quantile(W_values, c(0.025, 0.975), na.rm=T)

  # Extract the elements of the CI
  W.ci.lower <- round(CI[[1]],3)
  W.ci.upper <- round(CI[[2]],3)

  # create a vector to store the results to be later reported in the output table
  report.of.W <- paste0(W, " (95% CI ", W.ci.lower,"-", W.ci.upper, ")")

  # Use the cut() function to determine the effect size for W
  W.effect.size <- cut(W, breaks = breakpoints_Cohen, labels = labels_Cohen, right = TRUE)


  ##Cohen's w##
  w <- round(V * sqrt(min(nr, nc)-1),3)

  # Define breakpoints for effect sizes
  breakpoints_w <- c(-Inf, 0.10, 0.30, 0.50, Inf)
  labels_w <- c("(negligible effect)", "(small effect)", "(medium effect)", "(large effect)")

  # Use the cut() function to determine the effect size for w
  w_label <- cut(w, breaks = breakpoints_w, labels = labels_w, right = TRUE)


  ##Goodman-Kruskal's lambda##
  E1 <- n - max(sr)
  E2 <- sum(sc - apply(df, 2, max))
  lambda.row.dep <- round((E1-E2)/E1,3)

  E1 <- n - max(sc)
  E2 <- sum(sr - apply(df, 1, max))
  lambda.col.dep <- round((E1-E2)/E1,3)

  ##Goodman-Kruskal's lambda symmetric##
  max.col.wise <- apply(df, 2, max)
  max.row.wise <- apply(df, 1, max)
  max.row.tot <- max(sr)
  max.col.tot <- max(sc)
  lambda <- round(((sum(max.col.wise) + sum(max.row.wise)) - max.row.tot - max.col.tot) / ((2*n)-max.row.tot-max.col.tot),3)


  ## Define a function to calculate the corrected Lambda ##
  lambda_corrected <- function(table) {
    # Normalize the table if it contains counts
    if (max(table) > 1) {
      table <- table / sum(table)
    }

    # Calculating row and column sums
    p_iplus <- rowSums(table)
    p_plusj <- colSums(table)

    # Calculating maximums for rows and columns
    p_im <- apply(table, 1, max)
    p_mj <- apply(table, 2, max)

    # Calculating the asymmetric coefficients
    lambda_Y_given_X_K <- (sqrt(sum(p_im^2 / p_iplus)) - max(p_plusj)) / (1 - max(p_plusj))
    lambda_X_given_Y_K <- (sqrt(sum(p_mj^2 / p_plusj)) - max(p_iplus)) / (1 - max(p_iplus))

    # Calculating the symmetric coefficient
    M <- (sqrt(sum(p_im^2 / p_iplus)) + sqrt(sum(p_mj^2 / p_plusj)) - max(p_plusj) - max(p_iplus)) /
      (2 - max(p_plusj) - max(p_iplus))

    return(list(lambda_Y_given_X_K = lambda_Y_given_X_K,
                lambda_X_given_Y_K = lambda_X_given_Y_K,
                M = M))
  }

  # Put the above function to work and get the corrected Lambdas
  lambdas <- lambda_corrected(df)
  lambda.corrected.row.dep <- round(lambdas$lambda_Y_given_X_K, 3)
  lambda.corrected.col.dep <- round(lambdas$lambda_X_given_Y_K, 3)
  lambda.corrected.symm <- round(lambdas$M, 3)


  ##Goodman-Kruskal's tau##
  calc.tau <- function(x){
    tot.n.errors.rowwise <- sum(((sum(x)-rowSums(x))/sum(x))*rowSums(x))
    errors.col.wise <- matrix(nrow=nrow(x), ncol=ncol(x))
    for (i in 1:nrow(x)) {
      for (j in 1:ncol(x)) {
        errors.col.wise[i,j] <- ((colSums(x)[j] - x[i,j]) / colSums(x)[j]) * x[i,j]
      }
    }
    tau <- round((tot.n.errors.rowwise - sum(errors.col.wise)) / tot.n.errors.rowwise, 3)
    return(tau)
  }

  tau.row.dep <- calc.tau(t(df))
  tau.col.dep <- calc.tau(df)


  ##Odds ratio##
  # Check if the table is 2x2
  if (nr==2 & nc==2) {
    # Check for zeros along the diagonal of the 2x2 table
    if (df[1,1] * df[2,2] == 0 || df[1,2] * df[2, 1] == 0) {
      # In case of zeros along any of the diagonals, add 0.5 to every cell of the 2x2 table
      df.to.use <- df + 0.5
    } else {
      df.to.use <- df
    }
    or <- round(((df.to.use[1,1]*df.to.use[2,2]) / (df.to.use[1,2]*df.to.use[2,1])),3)
    se <- sqrt((1/df.to.use[1,1])+(1/df.to.use[1,2])+(1/df.to.use[2,1])+(1/df.to.use[2,2]))
    or.upper.ci <- round(exp(log(or)+1.96*se),3)
    or.lower.ci <- round(exp(log(or)-1.96*se),3)
    or.z <- (log(or)/se)
    or.p.value <- as.numeric(format((exp(-0.717*or.z-0.416*or.z*or.z)), scientific = T))

    or.p.to.report <- ifelse(or.p.value < 0.001, "< 0.001",
                             ifelse(or.p.value < 0.01, "< 0.01",
                                    ifelse(or.p.value < 0.05, "< 0.05",
                                           round(or.p.value,3))))

    report.of.or <- paste0(or, " (95% CI ", or.lower.ci,"-", or.upper.ci, "; p: ", or.p.to.report, ")")

  } else {
    report.of.or <- "-"
    or <- "-"
    or.lower.ci <- "-"
    or.upper.ci <- "-"
    or.p.value <- "-"
  }

  #Define a function to calculate the effect-size for the odds ratio
  #based on the thresholds suggested by Ferguson 2009.
  effect_size_from_OR <- function(OR) {

    # If the OR is less than 1, consider its reciprocal for categorization
    if (OR < 1) {
      OR <- 1 / OR
    }

    # Check the OR value
    if (OR < 2.0) {
      return("negligible")
    } else if (OR >= 2.0 && OR < 3.0) {
      return("small effect")
    } else if (OR >= 3.0 && OR < 4.0) {
      return("medium effect")
    } else {
      return("large effect")
    }
  }

  #Put the above function to work and work out the
  #odds ratio's effect size
  if (or == "-") {
    or.effect.size <- "-"
    or.effect.size.to.report <- ""
  } else {
    or.effect.size <- effect_size_from_OR(or)
    or.effect.size.to.report <- paste0(" (", or.effect.size,")")
  }

  ##Compute ORs for tables larger than 2x2
  # Check if the table is 2x2
  if (nr > 2 || nc > 2) {
    ORs.table <- round(calculate_odds_ratios(df),3)
    note.for.ORs <- "*Odds ratios are computed for adjacent rows and columns, representing the minimal set of odds ratios
    needed to describe all pairwise relationships in the contingency table.*"
  } else {
    ORs.table <- "-"
  }


  ##compute the power of the chi-square test##
  c <- qchisq(1-power.alpha, degrees.of.f)
  chi.sq.power <- 1 - pchisq(c, degrees.of.f, chisq.stat)
  chi.sq.power.report <- paste0(round(chi.sq.power,3), " (alpha: ", power.alpha,")")


  ##standardised residuals##
  stand.res <- round((df - exp.freq) / sqrt(exp.freq), 3)


  ##moment-corrected standardised residuals##
  mom.corr.stand.res <- round(stand.res / (sqrt((nr-1)*(nc-1)/(nr*nc))),3)

  #create a copy of the previous dataframe to be populated
  #with the values of the adjusted standardised residuals
  adj.stand.res <- stand.res


  ##adjusted standardised residuals##
  for (i in 1:nr) {
    for (j in 1:nc) {
      adj.stand.res[i,j] <- round(stand.res[i,j] / sqrt((1-sr[i]/n)*(1-sc[j]/n)), 3)
    }
  }


  ##Quetelet Index##
  QI <- round((df / exp.freq) - 1, 3)
  #define a note to be used later on as annotation for the tables of Quetelet indices
  note.for.QI <-  "*BLUE: noteworthy reduction in probability (< -0.50) <br> RED: noteworthy increase in probability (> 1.00)*"


  ##IJ association factor##
  IJ_factor <- round(df / exp.freq, 3)
  note.for.IJ <-  "*BLUE: noteworthy reduction in probability (< 0.50) <br> RED: noteworthy increase in probability (> 2.00)*"


  ##Adjusted Standardised Counts##
  #compute the standardised table using the ancillary function with custom marginals (all set to unity)
  standardisation.for.adj.counts <- standardize_table(df, marginal.type = NULL, custom.row.totals = rep(1, nr), custom.col.totals = rep(1 , nc))

  #subtract the grand mean from the standardised counts
  adj.stand.counts <- round(standardisation.for.adj.counts$table.stand - (sum(standardisation.for.adj.counts$table.stand)/(nr*nc)),3)

  # extract the iterations number to be used later on in the 'gt' table annotation
  adj.stand.counts.n.iterations <- standardisation.for.adj.counts$n.iterations

  #define a note to be used as annotation
  note.for.adj.stand.counts <-  paste0("*Standardisation converged at iteration number ", adj.stand.counts.n.iterations, ".*")


  ##Chi-squares-maximising Table##
  #use the ancillary function to get relevant data
  chi.sq.max.data <- maximize_chi_squared(df)

  #extract the chi-squares-maximising Table
  chi.sq.max.table <- chi.sq.max.data$max_table

  #extract the chi-squared value
  chi.sq.max.stat <- round(chi.sq.max.data$chi_squared_max,3)

  #compute the sqrt of the ratio between the observed chi2 and the max chi2
  #chisq.to.chisq.max <- round(sqrt(chisq.stat/chi.sq.max.stat),3)

  #compute Vmax
  V.max <- round(sqrt(chi.sq.max.stat / (n * min(nr-1, nc-1))), 3)

  #compute V corrected
  V.corr <- round(V/V.max,3)

  # Use the cut() function to determine the effect size
  V.corr.effect.size <- cut(V.corr, breaks = breakpoints_Cohen, labels = labels_Cohen, right = FALSE)

  #define a note to be used later as annotation in the rendered table
  note.for.chi.sq.max.table <-  paste0("*Chi-square value from the chi-square-maximising table: ", chi.sq.max.stat, ".*")


  ##cells' relative contribution in percent to the chi-sq##
  relative.contrib <- round((chisq.values / chisq.stat * 100),3)


  ##cells' absolute contribution in percent to the chi-sq##
  absolute.contrib <- round((chisq.values / n * 100),3)


  ##average cells' relative contribution to chi-sq statistic##
  #to be used as threshold to give different colours
  #in the table produced later on with gt
  average.relat.contr <- 100 / (nr*nc)

  #average cells' absolute contribution to chi-sq statistic
  #to be used as threshold to give different colours
  #in the table produced later on with gt
  average.abs.contr <- sum(absolute.contrib) / (nr*nc)


  ##plot of the permutation and Monte Carlo distribution of the chi-square statistic##
  if (graph==TRUE) {

    #conditionally set the layout in just one visualization
    if(oneplot==TRUE){
      m <- rbind(c(1,2))
      layout(m)
    }

    graphics::hist(permuted.chi.squared,
                   main="Permutation Chi-Squared Statistic Distribution",
                   sub=paste0("\nBalck dot: observed chi-squared statistic (", round(chisq.stat, 3), ")",
                              "\nDashed line: 95th percentile of the permutation chi-squared statistic distribution (", round(quantile(permuted.chi.squared, c(0.95)),3),
                              ")", "\np-value: ", round(p.value.permut,3), " (95% CI: ", round(lower_ci.perm,3), "-", round(upper_ci.perm,3), ";", " n. of simulated tables: ", B,")"),
                   xlab = "",
                   cex.main=0.80,
                   cex.sub=0.60)
    graphics::rug(chistat.mc, col = "#0000FF")
    graphics::points(x = chisq.stat, y=0, pch=20, col="black")
    graphics::abline(v = round(stats::quantile(chistat.mc, c(0.95)), 5), lty = 2, col = "blue")


    graphics::hist(chistat.mc,
                   main="Monte Carlo Chi-Squared Statistic Distribution",
                   sub=paste0("\nBalck dot: observed chi-squared statistic (", round(chisq.stat, 3), ")",
                              "\nDashed line: 95th percentile of the Monte Carlo chi-squared statistic distribution (", round(quantile(chistat.mc, c(0.95)),3),
                              ")", "\np-value: ", round(p.uppertail,3), " (95% CI: ", round(lower_ci,3), "-", round(upper_ci,3), ";", " n. of simulated tables: ", B,")"),
                   xlab = "",
                   cex.main=0.80,
                   cex.sub=0.60)
    graphics::rug(chistat.mc, col = "#0000FF")
    graphics::points(x = chisq.stat, y=0, pch=20, col="black")
    graphics::abline(v = round(stats::quantile(chistat.mc, c(0.95)), 5), lty = 2, col = "blue")

    #restore the original graphical device's settings if previously modified
    if(oneplot==TRUE){
      par(mfrow = c(1,1))
    }
  }


  ##put the internal 'suggest_chi_squared_method()' function to work##
  #and assign the output character vector to the decision_method object to be used
  #later on in the annotation of the first rendered table
  decision_method <- suggest_chi_squared_method(df)


  ##define a vector for the labels of the statistics to be returned##
  statistics <- c(
    "Expected Frequencies",
    "Smallest expected count",
    "Average expected count",
    "Hypothesis Tests",
    "Chi-Square Max",
    "Chi-Square Test",
    "Chi-Square Test (N-1)/N adjusted",
    "Chi-Square Test Permutation p-value",
    "Chi-Square Test Monte Carlo p-value",
    "G-Square Test",
    "Retrospective Power Analysis",
    "Chi-Square Test Power",
    "Margin-Free Association Measures",
    paste("Odds ratio", or.effect.size.to.report),
    paste("Yule's Q", Q.effect.size),
    paste("Yule's Y", Y.effect.size),
    "Chi-Square-Based Association Measures",
    paste("Phi signed ", phi_signed_effect_size),
    paste("Phi ", phi_effect_size),
    "Phi max ",
    paste("Phi corrected ", phi_corr_effect_size),
    paste("Cadj ", Cadj.effect.size),
    paste("Cramer's V ", V.effect.size),
    "Cramer's V's max",
    paste("Cramer's V corrected", V.corr.effect.size),
    paste("Cramer's V standardised ", V.stand.effect.size),
    "1-(Cramer's V / Cramer's V stand.)",
    paste("Cramer's V bias-corrected ", Vbc.effect.size),
    paste("Cohen's w ", w_label),
    paste("W ", W.effect.size),
    "PRE Association Measures",
    "Goodman-Kruskal's lambda (rows dependent)",
    "Goodman-Kruskal's lambda (columns dependent)",
    "Goodman-Kruskal's lambda (symmetric)",
    "Goodman-Kruskal's lambda corrected (rows dependent)",
    "Goodman-Kruskal's lambda corrected (columns dependent)",
    "Goodman-Kruskal's lambda corrected (symmetric)",
    "Goodmak-Kruskal's tau (rows dependent)",
    "Goodmak-Kruskal's tau (columns dependent)")


  ##define a vector for the statistics to be returned##
  values.to.report <- c(
    "",  # Empty text for the Expected Freq header
    round(min(exp.freq), 3),
    avrg.expt.count,
    "", # Hypoth tests: same as above
    chi.sq.max.stat,
    paste0(chisq.stat, " (df: ", degrees.of.f, "; p: ", p.to.report, ")"),
    paste0(chisq.adj, " (df: ", degrees.of.f, "; p: ", p.chisq.adj.to.report, ")"),
    paste0(round(p.value.permut,3), " (95% CI ", round(lower_ci.perm,3), "-", round(upper_ci.perm,3), ")"),
    paste0(round(p.uppertail,3), " (95% CI ", round(lower_ci,3), "-", round(upper_ci,3), ")"),
    paste0(Gsquared, " (p: ", p.Gsquared.to.report, ")"),
    "",  # Power: same as above
    chi.sq.power.report,
    "",  # Margin-free
    report.of.or,
    report.of.Q,
    report.of.Y,
    "",  # Chi-square-based
    phi.signed,
    phi,
    phi.max,
    phi.corr,
    Cadj,
    report.of.V,
    V.max,
    V.corr,
    V.stand,
    VtoV.stand,
    V.bc,
    w,
    report.of.W,
    "",  # PRE
    lambda.row.dep,
    lambda.col.dep,
    lambda,
    lambda.corrected.row.dep,
    lambda.corrected.col.dep,
    lambda.corrected.symm,
    tau.row.dep,
    tau.col.dep)


  ## Create the data frame with minimal, non-empty column names##
  val <- stat <- NULL
  report.df <- data.frame(stat = statistics, val = values.to.report)

  #define a vector of column names to be used later on
  #to conditionally give color to the cells text using 'gt'
  col.names.vect <- colnames(df)


  ##define the 'gt' elements for the output of the contingency table##
  #first, conditionally round the input table to eliminate the decimal places
  #that may have been been introduced in all the cells in an earlier step when replacing the 0s with 0.001
  #in order to be able to calculate G-squared
  if(any(df==0.001)){
    df <- round(df,0)
  }
  input.table.out <- as.data.frame(addmargins(as.matrix(df)))
  input.table.out <- gt::gt(input.table.out, rownames_to_stub = T)
  input.table.out <- gt::cols_align(input.table.out,align="center")
  input.table.out <- gt::tab_header(input.table.out,
                                    title = gt::md("**Analysis Report**"),
                                    subtitle = gt::md("*Observed frequencies*"))
  input.table.out <- gt::tab_source_note(input.table.out,
                                         md(paste0("Chi-Square Test: ", chisq.stat,
                                                   " (df: ", degrees.of.f, "; p: ", p.to.report,")
                                                   <br> Chi-Square Test (N-1)/N adj: ", chisq.adj, "(df: ", degrees.of.f, "; p: ", p.chisq.adj.to.report, ")
                                                   <br> Chi-Square Test Permutation p-value: ", round(p.value.permut,3), " (95% CI: ", round(lower_ci.perm,3), "-",round(upper_ci.perm,3), ")
                                                   <br> Chi-Square Test Monte Carlo p-value: ", round(p.uppertail,3), " (95% CI: ", round(lower_ci,3), "-",round(upper_ci,3), ")
                                                   <br> G-Square Test: ", Gsquared," (p: ", p.Gsquared.to.report, ")
                                                   <br><br> Chi-Square Test Power: ", chi.sq.power.report, "
                                                   <br><br> Yule's Q: ", report.of.Q,
                                                   "<br> Yule's Y: ", report.of.Y,
                                                   "<br> Odds Ratio: ", report.of.or,
                                                   "<br><br> Cramer's V: ", report.of.V,
                                                   "<br> Cramer's V corrected: ", V.corr,
                                                   "<br> Cramer's V standardised: ", V.stand,
                                                   "<br> 1-(Cramer's V / Cramer's V stand.): ", VtoV.stand,
                                                   "<br> W coefficient: ", report.of.W,
                                                   "<br><br> Goodman-Kruskal's lambda (rows depedent): ", lambda.row.dep,
                                                   "<br>Goodman-Kruskal's lambda (columns depedent): ", lambda.col.dep,
                                                   "<br><br> Number of cells with a significant standardised residual: ", sum(abs(stand.res) > z.crit.two.tailed),
                                                   "<br> Number of cells with a significant moment-corrected standardised residual: ", sum(abs(mom.corr.stand.res) > z.crit.two.tailed),
                                                   "<br> Number of cells with a significant adjusted standardised residual: ", sum(abs(adj.stand.res) > z.crit.two.tailed),
                                                   "<br> Significance of the residuals set at alpha ", round(alpha.adj,3), " (two-tailed z critical: ", round(z.crit.two.tailed,3), ")
                                                   <br><br>*", decision_method, ".*")))
  input.table.out <- gt::tab_options(input.table.out, source_notes.font.size=10, table.font.size=tfs, table.width = gt::pct(100))


  #define the 'gt' elements for the output table of the expected frequencies
  exp.freq.out <- gt::gt(as.data.frame(exp.freq), rownames_to_stub = T)
  exp.freq.out <- gt::cols_align(exp.freq.out,align="center")
  exp.freq.out <- gt::tab_header(exp.freq.out,
                                 title = gt::md("**Analysis Report**"),
                                 subtitle = gt::md("*Expected Frequencies*"))
  exp.freq.out <- gt::tab_options(exp.freq.out, table.font.size=tfs, table.width = gt::pct(100))


  #define the 'gt' elements for the output table of the chi-square values
  chisq.values.out <- gt::gt(as.data.frame(chisq.values), rownames_to_stub = T)
  chisq.values.out <- gt::cols_align(chisq.values.out,align="center")
  chisq.values.out <- gt::tab_header(chisq.values.out,
                                     title = gt::md("**Analysis Report**"),
                                     subtitle = gt::md("*Chi-Square Values*"))
  chisq.values.out <- gt::tab_options(chisq.values.out, table.font.size=tfs, table.width = gt::pct(100))


  #define the 'gt' elements for the output table of the relative contributions to the chi-square
  relative.contrib.out <- gt::gt(as.data.frame(relative.contrib), rownames_to_stub = T)
  relative.contrib.out <- gt::cols_align(relative.contrib.out,align="center")
  relative.contrib.out <- gt::tab_header(relative.contrib.out,
                                         title = gt::md("**Analysis Report**"),
                                         subtitle = gt::md("*Relative Contributions to the Chi-Square Statistic (in percent)*"))

  relative.contrib.out <- gt::tab_source_note(relative.contrib.out,
                                              md("*RED: larger-than-average contributions*"))
  relative.contrib.out <- gt::tab_options(relative.contrib.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))

  for(i in seq_along(col.names.vect)) {
    relative.contrib.out <- gt::tab_style(relative.contrib.out,
                                          style = gt::cell_text(color="red"),
                                          locations = gt::cells_body(
                                            columns = col.names.vect[i],
                                            rows = relative.contrib.out$`_data`[[col.names.vect[i]]] >= average.relat.contr))
  }


  #define the 'gt' elements for the output table of the absolute contributions to the chi-square
  absolute.contrib.out <- gt::gt(as.data.frame(absolute.contrib), rownames_to_stub = T)
  absolute.contrib.out <- gt::cols_align(absolute.contrib.out,align="center")
  absolute.contrib.out <- gt::tab_header(absolute.contrib.out,
                                         title = gt::md("**Analysis Report**"),
                                         subtitle = gt::md("*Absolute Contributions to the Shi-Square Statistic (in percent)*"))

  absolute.contrib.out <- gt::tab_source_note(absolute.contrib.out,
                                              md("*RED: larger-than-average contributions*"))
  absolute.contrib.out <- gt::tab_options(absolute.contrib.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))

  for(i in seq_along(col.names.vect)) {
    absolute.contrib.out <- gt::tab_style(absolute.contrib.out,
                                          style = gt::cell_text(color="red"),
                                          locations = gt::cells_body(
                                            columns = col.names.vect[i],
                                            rows = absolute.contrib.out$`_data`[[col.names.vect[i]]] >= average.abs.contr))
  }



  #define the 'gt' elements for the output table of the chi-sqaured-maximising cross-tab
  chi.sq.max.table.out <- as.data.frame(addmargins(as.matrix(chi.sq.max.table)))
  chi.sq.max.table.out <- gt::gt(chi.sq.max.table.out, rownames_to_stub = T)
  chi.sq.max.table.out <- gt::cols_align(chi.sq.max.table.out,align="center")
  chi.sq.max.table.out <- gt::tab_header(chi.sq.max.table.out,
                                         title = gt::md("**Analysis Report**"),
                                         subtitle = gt::md("*Chi-Square-Maximising Table*"))
  chi.sq.max.table.out <- gt::tab_source_note(chi.sq.max.table.out, md(note.for.chi.sq.max.table))
  chi.sq.max.table.out <- gt::tab_options(chi.sq.max.table.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))


  #define the 'gt' elements for the output table of the stand. residuals
  stand.res.out <- gt::gt(as.data.frame(stand.res), rownames_to_stub = T)
  stand.res.out <- gt::cols_align(stand.res.out,align="center")
  stand.res.out <- gt::tab_header(stand.res.out,
                                  title = gt::md("**Analysis Report**"),
                                  subtitle = gt::md("*Standardised Residuals*"))
  stand.res.out <- gt::tab_source_note(stand.res.out, md(note.for.residuals))
  stand.res.out <- gt::tab_options(stand.res.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))

  for(i in seq_along(col.names.vect)) {
    stand.res.out <- gt::tab_style(stand.res.out,
                                   style = gt::cell_text(color="red"),
                                   locations = gt::cells_body(
                                     columns = col.names.vect[i],
                                     rows = stand.res.out$`_data`[[col.names.vect[i]]] >= z.crit.two.tailed))
  }

  for(i in seq_along(col.names.vect)) {
    stand.res.out <- gt::tab_style(stand.res.out,
                                   style = gt::cell_text(color="blue"),
                                   locations = gt::cells_body(
                                     columns = col.names.vect[i],
                                     rows = stand.res.out$`_data`[[col.names.vect[i]]] <= (0-z.crit.two.tailed)))
  }


  #define the 'gt' elements for the output table of the moment-corrected stand. residuals
  mom.corr.stand.res.out <- gt::gt(as.data.frame(mom.corr.stand.res), rownames_to_stub = T)
  mom.corr.stand.res.out <- gt::cols_align(mom.corr.stand.res.out,align="center")
  mom.corr.stand.res.out <- gt::tab_header(mom.corr.stand.res.out,
                                           title = gt::md("**Analysis Report**"),
                                           subtitle = gt::md("*Moment-Corrected Standardised Residuals*"))
  mom.corr.stand.res.out <- gt::tab_source_note(mom.corr.stand.res.out, md(note.for.residuals))
  mom.corr.stand.res.out <- gt::tab_options(mom.corr.stand.res.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))

  for(i in seq_along(col.names.vect)) {
    mom.corr.stand.res.out <- gt::tab_style(mom.corr.stand.res.out,
                                            style = gt::cell_text(color="red"),
                                            locations = gt::cells_body(
                                              columns = col.names.vect[i],
                                              rows = mom.corr.stand.res.out$`_data`[[col.names.vect[i]]] >= z.crit.two.tailed))
  }

  for(i in seq_along(col.names.vect)) {
    mom.corr.stand.res.out <- gt::tab_style(mom.corr.stand.res.out,
                                            style = gt::cell_text(color="blue"),
                                            locations = gt::cells_body(
                                              columns = col.names.vect[i],
                                              rows = mom.corr.stand.res.out$`_data`[[col.names.vect[i]]] <= (0-z.crit.two.tailed)))
  }


  #define the 'gt' elements for the output table of the adjusted stand. residuals
  adj.stand.res.out <- gt::gt(as.data.frame(adj.stand.res), rownames_to_stub = T)
  adj.stand.res.out <- gt::cols_align(adj.stand.res.out,align="center")
  adj.stand.res.out <- gt::tab_header(adj.stand.res.out,
                                      title = gt::md("**Analysis Report**"),
                                      subtitle = gt::md("*Adjusted Standardised Residuals*"))
  adj.stand.res.out <- gt::tab_source_note(adj.stand.res.out, md(note.for.residuals))
  adj.stand.res.out <- gt::tab_options(adj.stand.res.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))

  for(i in seq_along(col.names.vect)) {
    adj.stand.res.out <- gt::tab_style(adj.stand.res.out,
                                       style = gt::cell_text(color="red"),
                                       locations = gt::cells_body(
                                         columns = col.names.vect[i],
                                         rows = adj.stand.res.out$`_data`[[col.names.vect[i]]] >= z.crit.two.tailed))
  }

  for(i in seq_along(col.names.vect)) {
    adj.stand.res.out <- gt::tab_style(adj.stand.res.out,
                                       style = gt::cell_text(color="blue"),
                                       locations = gt::cells_body(
                                         columns = col.names.vect[i],
                                         rows = adj.stand.res.out$`_data`[[col.names.vect[i]]] <= (0-z.crit.two.tailed)))
  }


  #define the 'gt' elements for the output table of the table of ORs
  if (nr > 2 || nc > 2) {
    ORs.out <- gt::gt(as.data.frame(ORs.table), rownames_to_stub = T)
    ORs.out <- gt::cols_align(ORs.out,align="center")
    ORs.out <- gt::tab_header(ORs.out,
                              title = gt::md("**Analysis Report**"),
                              subtitle = gt::md("*Independent Odds Ratios*"))
    ORs.out <- gt::tab_source_note(ORs.out, md(note.for.ORs))
    ORs.out <- gt::tab_options(ORs.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))
  }


  #define the 'gt' elements for the output table of the Quetelet Index
  QI.out <- gt::gt(as.data.frame(QI), rownames_to_stub = T)
  QI.out <- gt::cols_align(QI.out,align="center")
  QI.out <- gt::tab_header(QI.out,
                           title = gt::md("**Analysis Report**"),
                           subtitle = gt::md("*Quetelet Index*"))
  QI.out <- gt::tab_source_note(QI.out, md(note.for.QI))
  QI.out <- gt::tab_options(QI.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))

  for(i in seq_along(col.names.vect)) {
    QI.out <- gt::tab_style(QI.out,
                            style = gt::cell_text(color="red"),
                            locations = gt::cells_body(
                              columns = col.names.vect[i],
                              rows = QI.out$`_data`[[col.names.vect[i]]] > 1.00))
  }

  for(i in seq_along(col.names.vect)) {
    QI.out <- gt::tab_style(QI.out,
                            style = gt::cell_text(color="blue"),
                            locations = gt::cells_body(
                              columns = col.names.vect[i],
                              rows = QI.out$`_data`[[col.names.vect[i]]] < -0.50))
  }



  #define the 'gt' elements for the output table of the IJ association factor
  IJ.out <- gt::gt(as.data.frame(IJ_factor), rownames_to_stub = T)
  IJ.out <- gt::cols_align(IJ.out,align="center")
  IJ.out <- gt::tab_header(IJ.out,
                           title = gt::md("**Analysis Report**"),
                           subtitle = gt::md("*IJ Association Factor*"))
  IJ.out <- gt::tab_source_note(IJ.out, md(note.for.IJ))
  IJ.out <- gt::tab_options(IJ.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))

  for(i in seq_along(col.names.vect)) {
    IJ.out <- gt::tab_style(IJ.out,
                            style = gt::cell_text(color="red"),
                            locations = gt::cells_body(
                              columns = col.names.vect[i],
                              rows = IJ.out$`_data`[[col.names.vect[i]]] > 2.00))
  }

  for(i in seq_along(col.names.vect)) {
    IJ.out <- gt::tab_style(IJ.out,
                            style = gt::cell_text(color="blue"),
                            locations = gt::cells_body(
                              columns = col.names.vect[i],
                              rows = IJ.out$`_data`[[col.names.vect[i]]] < 0.50))
  }


  #define the 'gt' elements for the output table of the adjusted standardised counts
  adj.std.counts.out <- gt::gt(as.data.frame(adj.stand.counts), rownames_to_stub = T)
  adj.std.counts.out <- gt::cols_align(adj.std.counts.out,align="center")
  adj.std.counts.out <- gt::tab_header(adj.std.counts.out,
                                       title = gt::md("**Analysis Report**"),
                                       subtitle = gt::md("*Adjusted Standardised Counts*"))
  adj.std.counts.out <- gt::tab_source_note(adj.std.counts.out, md(note.for.adj.stand.counts))
  adj.std.counts.out <- gt::tab_options(adj.std.counts.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))


  #define the 'gt' elements for the output standardised table
  stand.table.out <- as.data.frame(addmargins(as.matrix(round(table.stand,3))))
  stand.table.out  <- gt::gt(stand.table.out, rownames_to_stub = T)
  stand.table.out  <- gt::cols_align(stand.table.out,align="center")
  stand.table.out  <- gt::tab_header(stand.table.out,
                                     title = gt::md("**Analysis Report**"),
                                     subtitle = gt::md("*Table Standardised via Iterative Proportional Fitting*"))

  #define a note to be used as annotation
  note.for.stand.table <-  paste0("*Standardisation converged at iteration number ", n.iterations, ".*")

  stand.table.out <- gt::tab_source_note(stand.table.out, md(note.for.stand.table))
  stand.table.out <- gt::tab_options(stand.table.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))


  #define the 'gt' elements for the output table of the analysis' results
  report.out <- gt::gt(report.df)

  # Hide column labels
  report.out <- gt::cols_label(report.out,
                               stat = "",
                               val = "")

  # Align values to the right
  report.out <- gt::cols_align(report.out, align = "right", columns = val)

  report.out <- gt::tab_header(report.out,
                               title = gt::md("**Analysis Report**"))

  # Style the section headers
  report.out <- gt::tab_style(
    report.out,
    style = list(
      gt::cell_text(weight = "bold", color = "black"),
      gt::cell_fill(color = "lightgrey")
    ),
    locations = gt::cells_body(
      columns = gt::everything(),
      rows = stat %in% c("Expected Frequencies","Hypothesis Tests", "Retrospective Power Analysis", "Margin-Free Association Measures", "Chi-Square-Based Association Measures", "PRE Association Measures")
    )
  )

  report.out <- gt::tab_source_note(report.out,
                                    md("Alternative thresholds for the strength of categorical association: <br>
                                           *Weak* association 0.00-0.10; *Moderate* association 0.11-0.30; *Strong* association 0.31-1.00
                                           <br> <br> after Healey J.F. 2013, *The Essentials of Statistics. A Tool for Social Research*, 3rd ed., Belmont: Wadsworth"))

  report.out <- gt::tab_options(report.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))


  #print the various 'gt' formatted tables defined in the above steps
  print(input.table.out)
  print(exp.freq.out)
  print(chisq.values.out)
  print(relative.contrib.out)
  print(absolute.contrib.out)
  print(chi.sq.max.table.out)
  print(stand.res.out)
  print(mom.corr.stand.res.out)
  print(adj.stand.res.out)
  if (nr > 2 || nc > 2) {
    print(ORs.out)
  }
  print(QI.out)
  print(IJ.out)
  print(adj.std.counts.out)
  print(stand.table.out)
  print(report.out)


  # Grouped results
  results <- list(
    input.table = list(
      "crosstab" = df
    ),
    chi.sq.maxim.table = list(
      "chi.sq.maximising.table" = chi.sq.max.table
    ),
    standardised.table = list(
      "standard.table" = round(table.stand,3)
    ),
    chi.sq.related.results = list(
      "exp.freq"=exp.freq,
      "smallest.exp.freq"=round(min(exp.freq),3),
      "avrg.exp.freq"=avrg.expt.count,
      "chisq.values"=chisq.values,
      "chisq.relat.contrib"=relative.contrib,
      "chisq.abs.contrib"=absolute.contrib,
      "chisq.statistic"=chisq.stat,
      "chisq.p.value"=p,
      "chisq.max"=chi.sq.max.stat,
      "chi.sq.power"=chi.sq.power,
      "chisq.adj"=chisq.adj,
      "chisq.adj.p.value"=p.chisq.adj,
      "chisq.p.value.perm"=p.value.permut,
      "chisq.p.value.perm CI lower boundary"=lower_ci.perm,
      "chisq.p.value.perm CI upper boundary"=upper_ci.perm,
      "chisq.p.value.MC"=p.uppertail,
      "chisq.p.value.MC CI lower boundary"=lower_ci,
      "chisq.p.value.MC CI upper boundary"=upper_ci
    ),
    G.square = list(
      "Gsq.statistic"=Gsquared,
      "Gsq.p.value"=p.Gsquared
    ),
    post.hoc = list(
      "stand.resid"=stand.res,
      "mom.corr.stand.resid"=mom.corr.stand.res,
      "adj.stand.resid"=adj.stand.res,
      "Quetelet.Index"=QI,
      "IJ.assoc.fact."=IJ_factor,
      "adj.stand.counts"=adj.stand.counts
    ),
    margin.free.assoc.measures = list(
      "Yule's Q"=Q,
      "Yule's Q CI lower boundary"=Q.lower,
      "Yule's Q CI upper boundary"=Q.upper,
      "Yule's Q p.value"=Q.p,
      "Yule's Y"=Y,
      "Yule's Y CI lower boundary"=Y.lower,
      "Yule's Y CI upper boundary"=Y.upper,
      "Yule's Y p.value"=Y.p,
      "Odds ratio"=or,
      "Odds ratio CI lower boundary"=or.lower.ci,
      "Odds ratio CI upper boundary"=or.upper.ci,
      "Odds ratio p.value"=or.p.value,
      "ORs"=ORs.table
    ),
    chi.sq.based.assoc.measures = list(
      "Phi.signed"=phi.signed,
      "Phi"=phi,
      "Phi.max"=phi.max,
      "Phi.corr"=phi.corr,
      "Cadj"=Cadj,
      "Cramer's V"=V,
      "Cramer's V CI lower boundary"=V.lower,
      "Cramer's V CI upper boundary"=V.upper,
      "Cramer's V max"=V.max,
      "Cramer's V corr"=V.corr,
      "Cramer's V standard."=V.stand,
      "1-(Cramer's V / Cramer's V stand.)"=VtoV.stand,
      "Cramer's Vbc"=V.bc,
      "w"=w,
      "W"=W,
      "W CI lower boundary"=W.ci.lower,
      "W CI upper boundary"=W.ci.upper
    ),
    PRE.assoc.measures = list(
      "lambda (rows dep.)"=lambda.row.dep,
      "lambda (cols dep.)"=lambda.col.dep,
      "lambda (symmetric)"=lambda,
      "lambda corrected (rows dep.)"=lambda.corrected.row.dep,
      "lambda corrected (cols dep.)"=lambda.corrected.col.dep,
      "lambda corrected (symmetric)"=lambda.corrected.symm,
      "tau (rows dep.)"=tau.row.dep,
      "tau (cols dep.)"=tau.col.dep)
  )


  # conditionally run the function to visualise the odds ratios
  if (plot.or==TRUE) {
    visualize_odds_ratios(data, reference.level = reference.level, row.level = row.level , or.alpha = or.alpha)
  }

  return(results)
}
