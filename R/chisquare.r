#' R function for chi-square and G-square test of independence, measures of association, and standardized/moment-corrected
#' standardized/adjusted standardized residuals
#'
#' The function performs the chi-square (and the G-square) test of independence on the input contingency table,
#' calculates various measures of categorical association, returns standardized, moment-corrected standardized, and adjusted
#' standardized residuals (with indication of their significance), and calculates relative and absolute contributions to the chi-square.
#' The p value associated to the chi-square statistic is also calculated on the basis of a permutation-based procedure.
#' Nicely-formatted output tables are rendered.
#'
#' The following \strong{measures of categorical associations} are produced by the function:
#'  \itemize{
#'   \item Phi (only for 2x2 tables)
#'   \item Phi signed (only for 2x2 tables)
#'   \item Yule's Q (and p value; only for 2x2 tables)
#'   \item Odds ratio (with 95perc confidence interval and p value; only for 2x2 tables)
#'   \item Adjusted contingency coefficient C
#'   \item Cramer's V (with 95perc confidence interval; indication of the magnitude of the effect size according to Cohen is provided for tables with up to 5 degrees of freedom)
#'   \item Bias-corrected Cramer's V (indication of the magnitude of the effect size according to Cohen is provided for tables with up to 5 degrees of freedom)
#'   \item Cohen's w (with indication of the magnitude of the effect size)
#'   \item Goodman-Kruskal's lambda (asymmetric)
#'   \item Goodman-Kruskal's lambda (symmetric)
#'   \item Goodman-Kruskal's tau (asymmetric)
#'   \item Goodman-Kruskal's gamma (and p value)
#'   \item Cohen's k (and 95perc confidence interval)
#' }
#'
#' The \strong{p value} of the observed chi-square statistic is also calculated on the basis of a \strong{permutation-based approach},
#' using \emph{B} random tables created under the Null Hypothesis of independence. For the rationale of this approach,
#' see for instance the description in Beh E.J., Lombardo R. 2014, Correspondence Analysis: Theory, Practice
#' and New Strategies, Chichester, Wiley: 62-64.\cr
#'
#' The \strong{permutation-based p value} is calculated as follows: \cr
#' \eqn{(1 + sum (chistat.perm > chisq.stat)) / (1 + B)}, where
#' \emph{chistat.perm} is a vector storing the B chi-square statistics generated under the Null Hypothesis, and
#' \emph{chisq.stat} is the observed chi-square statistic. For the logic of the calculation, see for example
#' Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016: 387.\cr
#'
#' The \strong{moment-corrected standardized residuals} are calculated as follows: \cr
#' \eqn{stand.res / (sqrt((nr-1)*(nc-1)/(nr*nc)))}, where \emph{stand.res} is each cell's standardized residual, \emph{nr} and
#' \emph{nc} are the number of rows and columns respectively; see
#' Garcia-Perez, MA, and Nunez-Anton, V (2003). Cellwise Residual Analysis in Two-Way Contingency Tables. Educational and Psychological Measurement, 63(5): 827.\cr
#'
#' The \strong{adjusted standardized residuals} are calculated as follows: \cr
#' \eqn{stand.res[i,j] / sqrt((1-sr[i]/n)*(1-sc[j]/n))}, where \emph{stand.res} is the standardized residual for cell \emph{ij},
#' \emph{sr} is the row sum for row \emph{i}, \emph{sc} is the column sum for column \emph{j}, and
#' \emph{n} is the table grand total. The \emph{adjusted standardized residuals} may prove useful since it has
#' been demonstrated that the standardized residuals tend to underestimate the significance of differences in small samples.
#' The adjusted standardized residuals correct that deficiency.\cr
#'
#' The \strong{significance of the residuals} (standardized, moment-corrected standardized, and adjusted standardized) is assessed using alpha 0.05 or, optionally
#' (by setting the parameter 'adj.alpha' to TRUE),
#' using an adjusted alpha calculated using the Sidak's method:\cr
#' \eqn{alpha.adj = 1-(1 - 0.05)^(1/(nr*nc))}, where \emph{nr} and \emph{nc} are the number of rows and columns in the table respectively. The adjusted
#' alpha is then converted into a critical two-tailed z value. See: Beasley TM and Schumacker RE (1995), Multiple Regression Approach to Analyzing Contingency
#' Tables: Post Hoc and Planned Comparison Procedures, The Journal of Experimental Education, 64(1): 86, 89.
#'
#' The \strong{cells' relative contribution (in percent) to the chi-square statistic} is calculated as:\cr
#' \eqn{chisq.values / chisq.stat * 100}, where \emph{chisq.values} and \emph{chisq.stat} are the chi-square
#' value in each individual cell of the table and the value of the chi-square statistic, respectively. The
#' \emph{average contribution} is calculated as \eqn{100 / (nr*nc)}, where \emph{nr} and \emph{nc} are the
#' number of rows and columns in the table respectively.\cr
#'
#' The \strong{cells' absolute contribution (in percent) to the chi-square statistic} is calculated as:\cr
#' \eqn{chisq.values / n * 100}, where \emph{chisq.values} and \emph{n} are the chi-square
#' value in each individual cell of the table and the table's grant total, respectively. The
#' \emph{average contribution} is calculated as sum of all the absolute contributions divided by the number of cells in
#' the table. For both the relative and absolute contributions to the chi-square, see
#' Beasley and Schumacker (1995): 90.\cr
#'
#' The calculation of the \strong{95perc confidence interval around Cramer's V} is based on Smithson M.J. (2003). Confidence Intervals,
#' Quantitative Applications in the Social Sciences Series, No. 140. Thousand Oaks, CA: Sage, 39-41, and builds on the R code made
#' available by the author on the web (http://www.michaelsmithson.online/stats/CIstuff/CI.html).\cr
#'
#' The \strong{bias-corrected Cramer's V} is based on Bergsma, W. (2013).
#' "A bias correction for Cramér's V and Tschuprow's T". Journal of the Korean Statistical Society. 42 (3): 323–328,
#' https://doi.org/10.1016%2Fj.jkss.2012.10.002\cr
#'
#' For the \strong{other measures of categorical association} provided by the function, see for example
#' Sheskin, D. J. (2011). Handbook of Parametric and Nonparametric Statistical Procedures, Fifth Edition (5th ed.). Chapman and Hall/CRC: 675-679,
#' 1415-1427.\cr
#'
#' \strong{Note} that:\cr
#' -the \strong{Phi} coefficient is based on the chi-square statistic as per Sheskin's equation 16.21, whereas the
#' \strong{Phi signed} is after Sheskin's equation 16.20;\cr
#'
#' -the \strong{2-sided p value of Yule's Q} is calculated following Sheskin's equation 16.24;\cr
#'
#' -\strong{Cohen's w} is calculate as \eqn{V * sqrt(min(nr, nc)-1)}, where \emph{V} is Cramer's V, and  \emph{nr} and  \emph{nc}
#' are the number of rows and columns respectively; see Sheskin 2011: 679;\cr
#'
#' -in the output table reporting the result of the chi-square test and the various measures of association, the
#' magnitude of the \strong{effect size} according to \emph{Cohen}'s guidelines is reported for Cramer's V, Cramer's V biase-corrected, and for Cohen's w;
#' the effect size for the former two coefficients is only reported for tables featuring up to 5 degrees of freedom
#' (see Cohen, J. (1988). Statistical power analysis for the behavioral sciences (2nd ed). Hillsdale, N.J: L. Erlbaum Associates);\cr
#'
#' -the \strong{2-tailed p value} of \strong{Goodman-Kruskal's gamma} is based on the
#' associated z-score calculated as per Sheskin's equation 32.2;\cr
#'
#' -the \strong{symmetric} version of \strong{Goodman-Kruskal's lambda} is calculated
#' as per Reynolds, H. T. (1984). Analysis of Nominal Data (Quantitative Applications in the Social Sciences) (1st ed.). SAGE Publications: 55-57;\cr
#'
#' -\strong{Goodman-Kruskal's tau} is calculated as per Reynolds 1984: 57-60;\cr
#'
#' -\strong{Cohen's k} is calculated as per Sheskin 2011: 688-689 (equation 16.30).\cr
#'
#' About the \strong{plot of odds ratios} for 2xk table, where k>2:\cr
#' by setting the \code{plor.or} parameter to \code{TRUE}, a plot showing the odds ratios and their 95percent confidence interval will be rendered.
#' The confidence level can be modified via the \code{or.alpha} parameter. The odds ratios are calculated for the column levels, and one of them
#' is to be selected by the user as a reference for comparison via the \code{reference_level} parameter (set to 1 by default).
#' Also, the user may want to select the row category to which the calculation of the odds ratios makes reference (using the \code{row_category} parameter,
#' which is set to 1 by default).\cr
#'
#' To better understand the rationale, consider the following example, which rests on the famous Titanic data:\cr
#'
#' Create a 2x3 contingency table:\cr
#' mytable <- matrix(c(123, 158, 528, 200, 119, 181), nrow = 2, byrow = TRUE)\cr
#' colnames(mytable) <- c("1st", "2nd", "3rd")\cr
#' rownames(mytable) <- c("Died", "Survived")\cr
#'
#' Now, we perform the test and visualise the odds ratios:\cr
#' chisquare(mytable, plot.or=TRUE, reference_level=1, row_category=1)\cr
#'
#' In the rendered plot, we can see the odds ratios and confidence intervals for the second and third column level
#' (i.e., 2nd class and 3rd class) because the first column level has been selected as reference level. The odds ratios are calculated
#' making reference to the first row category (i.e., \emph{Died}). From the plot, we can see that, compared to the 1st class,
#' passenger on the 2nd class have 2.16 times larger odds of dying; passengers on the 3rd class have 4.74 times larger odds of dying
#' compared to the 1st class.\cr
#'
#' Note that if we set the \code{row_category} parameter to \code{2}, we make reference to the second row category, i.e. \emph{Survived}:\cr
#' chisquare(mytable, plot.or=TRUE, reference_level=1, row_category=2)\cr
#'
#' In the plot, we can see that passengers in the 2nd class have 0.46 times the odds of surviving of passengers in the 1st class, while
#' passengers from the 3rd class have 0.21 times the odds of surviving of those who were travelling in the 1st class.\cr
#'
#'
#' @param data Dataframe containing the input contingency table.
#' @param B  Number of simulated tables to be used to calculate the Monte Carlo-based p value (999 by default).
#' @param plot.or Takes TRUE or FALSE (default) if the user wants a plot of the odds ratios to be rendered (only for 2xk tables, where k>2).
#' @param reference_level The index of the column reference level for odds ratio calculations (default: 1).
#' The user must select the column level to serve as the reference level (only for 2xk tables, where k>2).
#' @param row_category The index of the row category to be used in odds ratio calculations (1 or 2; default: 1).
#' The user must select the row level to which the calculation of the odds ratios make reference (only for 2xk tables, where k>2).
#' @param or.alpha The significance level used for the odds ratios' confidence intervals (default: 0.05).
#' @param adj.alpha  Takes TRUE or FALSE (default) if the user wants or does not want the significance level of the
#' residuals (both standarized and adjusted standardized) to be corrected using the Sidak's adjustment method (see Details).
#' @param format Takes \emph{short} (default) if the dataset is a dataframe storing a contingency table; if the
#' input dataset is a dataframe storing two columns that list the levels of the two categorical variables,
#' \emph{long} will preliminarily cross-tabulate the levels of the categorical variable in the 1st column against
#' the levels of the variable stored in the 2nd column.
#' @param graph Takes TRUE or FALSE (default) if the user wants or does not want to chart the distribution of the permuted
#'  chi-square statistic accross the number of randomized tables set by the B parameter.
#' @param tfs Numerical value to set the size of the font used in the main body of the various output tables (14 by default).
#'
#'
#' @return The function produces \strong{optional charts} (distribution of the permuted chi-square statistic
#' and a plot of the odds ratios between a reference column level and the other ones, the latter only for 2xk tables), and
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
#'   \item Standardized residuals (RED for large significant residuals, BLUE for small significant residuals)
#'   \item Moment-corrected standardized residuals (colour same as above)
#'   \item Adjusted standardized residuals (colour same as above)
#'   \item Table of output statistics, p values, and association measures
#' }
#'
#' Also, the function returns a \strong{list containing the following elements}:
#'  \itemize{
#'   \item \emph{crosstab}: input contingency table
#'   \item \emph{exp.freq.}: table of expected frequencies
#'   \item \emph{chisq.values}: cells' chi-square value
#'   \item \emph{chisq.relat.contrib.}: cells' relative contribution (in percent) to the chi-square statistic
#'   \item \emph{chisq.abs.contrib.}: cells' absolute contribution (in percent) to the chi-square statistic
#'   \item \emph{chisq.statistic}: observed chi-square value
#'   \item \emph{chisq.p.value}: p value of the chi-square statistic
#'   \item \emph{chisq.p.value.perm.}: p value based on B permuted tables
#'   \item \emph{Gsq.statistic}: observed G-square value
#'   \item \emph{Gsq.p.value}: p value of the G-square statistic
#'   \item \emph{stand.resid.}: table of standardized residuals
#'   \item \emph{mom.corr.stand.resid.}: table of moment-corrected standardized residuals
#'   \item \emph{adj.stand.resid.}: table of adjusted standardized residuals
#'   \item \emph{Phi}: Phi coefficient (only for 2x2 tables)
#'   \item \emph{Phi signed}: signed Phi coefficient (only for 2x2 tables)
#'   \item \emph{Yule's Q}: Q coefficient (only for 2x2 tables)
#'   \item \emph{Yule's Q p.value}: 2-tailed p value of Yule's Q
#'   \item \emph{Odds ratio}: odds ratio (only for 2x2 tables)
#'   \item \emph{Odds ratio CI lower boundary}: lower boundary of the 95perc CI
#'   \item \emph{Odds ratio CI upper boundary}: upper boundary of the 95perc CI
#'   \item \emph{Odds ratio p.value}: p value of the odds ratio
#'   \item \emph{Cadj}: adjusted contigency coefficient C
#'   \item \emph{Cramer's V}: Cramer's V coefficient
#'   \item \emph{Cramer's V CI lower boundary}: lower boundary of the 95perc CI
#'   \item \emph{Cramer's V CI upper boundary}: upper boundary of the 95perc CI
#'   \item \emph{Cramer's Vbc}: bias-corrected Cramer's V coefficient
#'   \item \emph{w}: Cohen's w
#'   \item \emph{lambda (rows dep.)}: Goodman-Kruskal's lambda coefficient (considering the rows being the dependent variable)
#'   \item \emph{lambda (cols dep.)}: Goodman-Kruskal's lambda coefficient (considering the columns being the dependent variable)
#'   \item \emph{lambda.symmetric}: Goodman-Kruskal's symmetric lambda coefficient
#'   \item \emph{tau (rows dep.)}: Goodman-Kruskal's tau coefficient (considering the rows being the dependent variable)
#'   \item \emph{tau (cols dep.)}: Goodman-Kruskal's tau coefficient (considering the columns being the dependent variable)
#'   \item \emph{gamma}: Goodman-Kruskal's gamma coefficient
#'   \item \emph{gamma.p.value}: 2-sided p value for the Goodman-Kruskal's gamma coefficient
#'   \item \emph{k}: Cohen'k
#'   \item \emph{k CI lower boundary}: lower boundary of the 95perc CI
#'   \item \emph{k CI upper boundary}: upper boundary of the 95perc CI
#' }
#'
#' \strong{Note} that the \emph{p values} returned in the above list are expressed in scientific notation, whereas the ones reported in the
#' output table featuring the tests' result and measures of association are reported as broken down into classes (e.g., <0.05, or <0.01, etc).\cr
#'
#' The \strong{following examples}, which use in-built datasets, can be run to fimiliarise with the function:\cr
#'
#' -perform the test on the in-built 'social_class' dataset\cr
#' result <- chisquare(social_class)\cr
#'
#' -perform the test on a 2x2 subset of the 'diseases' dataset\cr
#' mytable <- diseases[3:4,1:2]\cr
#' result <- chisquare(mytable)\cr
#'
#' -perform the test on a 2x2 subset of the 'safety' dataset\cr
#' mytable <- safety[c(4,1),c(1,6)]\cr
#' result <- chisquare(mytable)\cr
#'
#' -build a toy dataset in 'long' format (gender vs. opinion about death sentence)\cr
#' mytable <- data.frame(GENDER=c(rep("F", 360), rep("M", 340)),\cr
#' OPINION=c(rep("oppose", 235),\cr
#'          rep("favour", 125),\cr
#'          rep("oppose", 160),\cr
#'          rep("favour", 180)))\cr
#'
#' -perform the test specifying that the input table is in 'long' format\cr
#' result <- chisquare(mytable, format="long")\cr
#'
#'
#' @keywords chiperm
#'
#' @export
#'
#' @importFrom gt gt cols_align tab_header md tab_source_note tab_options pct
#' @importFrom stats pchisq addmargins r2dtable pnorm quantile qnorm
#' @importFrom graphics abline points hist rug
#'
#'
#' @examples
#' # Perform the test on the in-built 'social_class' dataset
#' result <- chisquare(social_class, B=99)
#'
#'
#'
chisquare <- function(data, B=999, plot.or= FALSE, reference_level = 1, row_category = 1, or.alpha = 0.05, adj.alpha=FALSE, format="short", graph=FALSE, tfs=14){
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

  #if alpha.adj is TRUE, adjust the value of alpha using the Sidak's method,
  #and calculate the two-tailed critical value to be used as threshold for the
  #significance of the standardized and adjusted standardized residuals
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
  #used later on to calculate the permuted chi-sq statistic
  calc <- function(x){
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
  #the chi-sq values, and the expected frequencies to be used later on
  chisq.stat <- calc(df)$chisq.stat
  chisq.values <- round(calc(df)$chisq.values,3)
  exp.freq <- round(calc(df)$exp.freq,3)

  #p value
  p <- as.numeric(format(pchisq(q=chisq.stat, df=degrees.of.f, lower.tail=FALSE), scientific = T))

  p.to.report <- ifelse(p < 0.001, "< 0.001",
                        ifelse(p < 0.01, "< 0.01",
                               ifelse(p < 0.05, "< 0.05",
                                      round(p, 3))))

  #calculate the chi-sq statistic B times, using the above-defined function
  extract <- function(x)  calc(x)$chisq.stat
  chistat.perm <- sapply(r2dtable(B, sr, sc), extract)

  #calculate the p value of the observed chi-sq on the basis of the
  #B permuted chi-sq statistics
  p.uppertail <- (1 + sum (chistat.perm > chisq.stat)) / (1 + B)

  p.uppertail.to.report <- ifelse(p.uppertail < 0.001, "< 0.001",
                        ifelse(p.uppertail < 0.01, "< 0.01",
                               ifelse(p.uppertail < 0.05, "< 0.05",
                                      round(p.uppertail, 3))))

  #G-squared
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

  #Phi
  if (nr==2 & nc==2) {
    phi <- round(sqrt(chisq.stat / n), 3)
    a <- as.numeric(df[1,1])
    b <- as.numeric(df[1,2])
    c <- as.numeric(df[2,1])
    d <- as.numeric(df[2,2])
    phi.signed <- round(((a*d)-(b*c)) / sqrt((a+b)*(c+d)*(a+c)*(b+d)),3)
  } else {
    phi <- "-"
    phi.signed <- "-"
  }

  #Yule's Q
  if (nr==2 & nc==2) {
    Q <- round(((df[1,1]*df[2,2]) - (df[1,2]*df[2,1])) /  ((df[1,1]*df[2,2]) + (df[1,2]*df[2,1])),3)
    Q.z <- Q / sqrt((1/4)*(1-Q^2)^2*(1/df[1,1]+1/df[1,2]+1/df[2,1]+1/df[2,2]))
    Q.p <- as.numeric(format(2*pnorm(q=abs(Q.z), lower.tail=FALSE), scientific=T))
    Q.p.to.report <- ifelse(Q.p < 0.001, "< 0.001",
                            ifelse(Q.p < 0.01, "< 0.01",
                                   ifelse(Q.p < 0.05, "< 0.05",
                                          round(Q.p, 3))))

    report.of.Q <- paste0(Q, " (p: ", Q.p.to.report, ")")
  } else {
    Q <- "-"
    Q.p <- "-"
    report.of.Q <- "-"
  }

  #adjusted contingency coeffcient
  C <- sqrt(chisq.stat / (n + chisq.stat))
  Cmax <- sqrt((min(nr,nc)-1) / min(nr,nc))
  Cadj <- round(C / Cmax,3)

  #Cramer's V
  V <- round(sqrt(chisq.stat / (n * min(nr-1, nc-1))), 3)

  #effect size for Cramer's V
  if (degrees.of.f <= 5 ) {
    if (degrees.of.f == 1){
      V.effect.size <- ifelse(V<=0.10, "negligible",
                              ifelse(V <= 0.30, "small effect",
                                     ifelse(V <= 0.50, "medium effect", "large effect")))
    }
    if (degrees.of.f == 2){
      V.effect.size <- ifelse(V<=0.07, "negligible",
                              ifelse(V<= 0.21, "small effect",
                                     ifelse(V <= 0.35, "medium effect", "large effect")))
    }
    if (degrees.of.f == 3){
      V.effect.size <- ifelse(V<=0.06, "negligible",
                              ifelse(V<= 0.17, "small effect",
                                     ifelse(V <= 0.29, "medium effect", "large effect")))
    }
    if (degrees.of.f == 4){
      V.effect.size <- ifelse(V<=0.05, "negligible",
                              ifelse(V<= 0.15, "small effect",
                                     ifelse(V <= 0.25, "medium effect", "large effect")))
    }
    if (degrees.of.f == 5){
      V.effect.size <- ifelse(V<=0.05, "negligible",
                              ifelse(V<= 0.13, "small effect",
                                     ifelse(V <= 0.22, "medium effect", "large effect")))
    }
    V.effect.size.to.report <- paste0(" (", V.effect.size,")")
  } else{
    V.effect.size.to.report <- ""
  }

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


  #bias-corrected Cramer's V
  phi.new <- max(0, (chisq.stat / n) - ((nc-1)*(nr-1)/(n-1)))
  k.new  <-  nc-((nc-1)^2 / (n-1))
  r.new  <-  nr-((nr-1)^2 / (n-1))
  V.bc <- round(sqrt(phi.new / min(k.new-1, r.new-1)), 3)

  #effect size for bias-corrected Cramer's V
  if (degrees.of.f <= 5 ) {

    if (degrees.of.f == 1){
      Vbc.effect.size <- ifelse(V.bc <= 0.10, "negligible",
                              ifelse(V.bc <= 0.30, "small effect",
                                     ifelse(V.bc <= 0.50, "medium effect", "large effect")))
    }

    if (degrees.of.f == 2){
      Vbc.effect.size <- ifelse(V.bc <= 0.07, "negligible",
                              ifelse(V.bc <= 0.21, "small effect",
                                     ifelse(V.bc <= 0.35, "medium effect", "large effect")))
    }

    if (degrees.of.f == 3){
      Vbc.effect.size <- ifelse(V.bc <= 0.06, "negligible",
                              ifelse(V.bc <= 0.17, "small effect",
                                     ifelse(V.bc <= 0.29, "medium effect", "large effect")))
    }
    if (degrees.of.f == 4){
      Vbc.effect.size <- ifelse(V.bc <= 0.05, "negligible",
                              ifelse(V.bc <= 0.15, "small effect",
                                     ifelse(V.bc <= 0.25, "medium effect", "large effect")))
    }
    if (degrees.of.f == 5){
      Vbc.effect.size <- ifelse(V.bc <= 0.05, "negligible",
                              ifelse(V.bc <= 0.13, "small effect",
                                     ifelse(V.bc <= 0.22, "medium effect", "large effect")))
    }
    Vbc.effect.size.to.report <- paste0(" (", Vbc.effect.size,")")

  } else{
    Vbc.effect.size.to.report <- ""
  }


  #Cohen's w
  w <- round(V * sqrt(min(nr, nc)-1),3)
  w_label <- ifelse(w <= 0.10, "negligible",
                    ifelse(w <= 0.30, "small effect",
                            ifelse(w <= 0.50, "medium effect", "large effect")))

  #Goodman-Kruskal's lambda
  E1 <- n - max(sr)
  E2 <- sum(sc - apply(df, 2, max))
  lambda.row.dep <- round((E1-E2)/E1,3)

  E1 <- n - max(sc)
  E2 <- sum(sr - apply(df, 1, max))
  lambda.col.dep <- round((E1-E2)/E1,3)

  #Goodman-Kruskal's lambda symmetric
  max.col.wise <- apply(df, 2, max)
  max.row.wise <- apply(df, 1, max)
  max.row.tot <- max(sr)
  max.col.tot <- max(sc)
  lambda <- round(((sum(max.col.wise) + sum(max.row.wise)) - max.row.tot - max.col.tot) / ((2*n)-max.row.tot-max.col.tot),3)

  #Goodman-Kruskal's gamma
  #define the function to calculate the coefficient
  calc.gamma <- function(x){
    n <- nrow(x)
    m <- ncol(x)
    pi.c <- pi.d <- matrix(0, nrow = n, ncol = m)
    row.x <- row(x)
    col.x <- col(x)
    for (i in 1:n) {
      for (j in 1:m) {
        pi.c[i, j] <- sum(x[row.x < i & col.x < j]) + sum(x[row.x > i & col.x > j])
        pi.d[i, j] <- sum(x[row.x < i & col.x > j]) + sum(x[row.x > i & col.x < j])
      }
    }
    C <- sum(pi.c * x)/2
    D <- sum(pi.d * x)/2
    GK.gamma <- (C - D)/(C + D)
    GK.gamma.z <- GK.gamma * sqrt((C+D)/(sum(x)*(1-GK.gamma^2)))
    #2-tailed p value from z
    GK.gamma.p <- 2*pnorm(q=abs(GK.gamma.z), lower.tail=FALSE)
    return(list("GK.gamma"=GK.gamma,
                "GK.gamma.p"=GK.gamma.p))
  }
  #put the above function to work to calculate gamma and its p value
  gamma.coeff <- round(calc.gamma(df)$GK.gamma,3)
  gamma.p <- as.numeric(format(calc.gamma(df)$GK.gamma.p, scientific = T))

  gamma.p.to.report <- ifelse(gamma.p < 0.001, "< 0.001",
                        ifelse(gamma.p < 0.01, "< 0.01",
                               ifelse(gamma.p < 0.05, "< 0.05",
                                      round(gamma.p,3))))

  #Goodman-Kruskal's tau
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

  #Odds ratio
  if (nr==2 & nc==2) {
    or <- round(((df[1,1]*df[2,2]) / (df[1,2]*df[2,1])),3)
    se <- sqrt((1/df[1,1])+(1/df[1,2])+(1/df[2,1])+(1/df[2,2]))
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

  #Cohen's k
  sum.observed.diag <- sum(diag(as.matrix(df)))
  sum.expected.diag <- sum(diag(as.matrix(exp.freq)))
  cohen.k <- round((sum.observed.diag - sum.expected.diag) / (n - sum.expected.diag),3)
  SD.k <- round(sqrt((sum.observed.diag * (n-sum.observed.diag)) / (n * (n-sum.expected.diag)^2)),3)
  k.lower.ci <- round(cohen.k - (1.96 * SD.k),3)
  k.upper.ci <- round(cohen.k + (1.96 * SD.k),3)
  report.of.k <- paste0(cohen.k, " (95% CI ", k.lower.ci,"-", k.upper.ci, ")")

  #standardized residuals
  stand.res <- round((df - exp.freq) / sqrt(exp.freq), 3)

  #moment-corrected standardized residuals
  mom.corr.stand.res <- round(stand.res / (sqrt((nr-1)*(nc-1)/(nr*nc))),3)

  #create a copy of the previous dataframe to be populated
  #with the values of the adjusted standardized residuals
  adj.stand.res <- stand.res

  #adjusted standardized residuals
  for (i in 1:nr) {
    for (j in 1:nc) {
      adj.stand.res[i,j] <- round(stand.res[i,j] / sqrt((1-sr[i]/n)*(1-sc[j]/n)), 3)
    }
  }

  #cells' relative contribution in percent to the chi-sq
  relative.contrib <- round((chisq.values / chisq.stat * 100),3)

  #cells' absolute contribution in percent to the chi-sq
  absolute.contrib <- round((chisq.values / n * 100),3)

  #average cells' relative contribution to chi-sq statistic
  #to be used as threshold to give different colours
  #in the table produced later on with gt
  average.relat.contr <- 100 / (nr*nc)

  #average cells' absolute contribution to chi-sq statistic
  #to be used as threshold to give different colours
  #in the table produced later on with gt
  average.abs.contr <- sum(absolute.contrib) / (nr*nc)

  #plot of the permuted distribution of the chi-square statistic
  if (graph==TRUE) {
    graphics::hist(chistat.perm,
                   main="Permuted Distribution of the Chi-square Statistic",
                   sub=paste0("\nBalck dot: observed chi-square statistic (", round(chisq.stat, 3), ")",
                              "\nDashed line: 95th percentile of the permuted chi-square statistic distribution (", round(quantile(chistat.perm, c(0.95)),3),
                              ")", "\np value: ", p.uppertail.to.report, " (n. of permutations: ", B,")"),
                   xlab = "",
                   cex.main=0.85,
                   cex.sub=0.70)
    graphics::rug(chistat.perm, col = "#0000FF")
    graphics::points(x = chisq.stat, y=0, pch=20, col="black")
    graphics::abline(v = round(stats::quantile(chistat.perm, c(0.95)), 5), lty = 2, col = "blue")
  }

  #define a vector for the labels of the statistics to be returned
  statistics <- c("Chi-square",
                  "G-square",
                  "Phi",
                  "Phi signed",
                  "Yule's Q",
                  "Odds ratio",
                  "Cadj",
                   paste0("Cramer's V", V.effect.size.to.report),
                   paste0("Cramer's V bias-corrected", Vbc.effect.size.to.report),
                   paste0("Cohen's w", " (", w_label, ")"),
                  "Goodman-Kruskal's lambda (rows dependent)",
                  "Goodman-Kruskal's lambda (columns dependent)",
                  "Goodman-Kruskal's lambda (symmetric)",
                  "Goodmak-Kruskal's tau (rows dependent)",
                  "Goodmak-Kruskal's tau (columns dependent)",
                  "Goodmak-Kruskal's gamma",
                  "Cohen's k")

  #define a vector for the statistics to be returned
  values.to.report <-c(paste0(chisq.stat, " (df: ", degrees.of.f, "; p: ", p.to.report, "; perm. p: ", p.uppertail.to.report, ")"),
                       paste0(Gsquared, " (p: ", p.Gsquared.to.report, ")"),
                       phi,
                       phi.signed,
                       report.of.Q,
                       report.of.or,
                       Cadj,
                       report.of.V,
                       V.bc,
                       w,
                       lambda.row.dep,
                       lambda.col.dep,
                       lambda,
                       tau.row.dep,
                       tau.col.dep,
                       paste0(gamma.coeff, " (p: ", gamma.p.to.report, ")"),
                       report.of.k)

  #create a dataframe storing the two above-defined vectors
  #this will be formatted later on according to the 'gt' settings
  report.df <- data.frame(STATISTIC=statistics, VALUE=values.to.report)


  #define a vector of column names to be used later on
  #to conditionally give color to the cells text using 'gt'
  col.names.vect <- colnames(df)

  #define the 'gt' elements for the output of the contingency table
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
                                 title = gt::md("**Analysis report**"),
                                 subtitle = gt::md("*Observed frequencies*"))
  input.table.out <- gt::tab_source_note(input.table.out,
                                    md(paste0("*Chi-square: ", chisq.stat,
                                    " (df: ", degrees.of.f, "; p: ", p.to.report,") <br>
                                    G-square: ", Gsquared," (p: ", p.Gsquared.to.report, ") <br>
                                    Cramer's V: ", report.of.V,
                                    "<br> Goodman-Kruskal's lambda (symmetric): ", lambda,
                                    "<br> Goodman-Kruskal's gamma: ", gamma.coeff,
                                    "<br> Number of cells with a significant standardized residual: ", sum(abs(stand.res) > z.crit.two.tailed),
                                    "<br> Number of cells with a significant moment-corrected standardized residual: ", sum(abs(mom.corr.stand.res) > z.crit.two.tailed),
                                    "<br> Number of cells with a significant adjusted standardized residual: ", sum(abs(adj.stand.res) > z.crit.two.tailed),
                                    "<br> Significance of the residuals set at alpha ", round(alpha.adj,3), " (two-tailed z critical: ", round(z.crit.two.tailed,3), ")*")))
  input.table.out <- gt::tab_options(input.table.out, source_notes.font.size=10, table.font.size=tfs, table.width = gt::pct(100))


  #define the 'gt' elements for the output table of the expected frequencies
  exp.freq.out <- gt::gt(as.data.frame(exp.freq), rownames_to_stub = T)
  exp.freq.out <- gt::cols_align(exp.freq.out,align="center")
  exp.freq.out <- gt::tab_header(exp.freq.out,
                             title = gt::md("**Analysis report**"),
                             subtitle = gt::md("*Expected frequencies*"))
  exp.freq.out <- gt::tab_options(exp.freq.out, table.font.size=tfs, table.width = gt::pct(100))


  #define the 'gt' elements for the output table of the chi-square values
  chisq.values.out <- gt::gt(as.data.frame(chisq.values), rownames_to_stub = T)
  chisq.values.out <- gt::cols_align(chisq.values.out,align="center")
  chisq.values.out <- gt::tab_header(chisq.values.out,
                                 title = gt::md("**Analysis report**"),
                                 subtitle = gt::md("*Chi-square values*"))
  chisq.values.out <- gt::tab_options(chisq.values.out, table.font.size=tfs, table.width = gt::pct(100))


  #define the 'gt' elements for the output table of the relative contributions to the chi-square
  relative.contrib.out <- gt::gt(as.data.frame(relative.contrib), rownames_to_stub = T)
  relative.contrib.out <- gt::cols_align(relative.contrib.out,align="center")
  relative.contrib.out <- gt::tab_header(relative.contrib.out,
                                     title = gt::md("**Analysis report**"),
                                     subtitle = gt::md("*Relative contributions to the chi-square statistic (in percent)*"))

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
                                         title = gt::md("**Analysis report**"),
                                         subtitle = gt::md("*Absolute contributions to the chi-square statistic (in percent)*"))

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


  #define the 'gt' elements for the output table of the stand. residuals
  stand.res.out <- gt::gt(as.data.frame(stand.res), rownames_to_stub = T)
  stand.res.out <- gt::cols_align(stand.res.out,align="center")
  stand.res.out <- gt::tab_header(stand.res.out,
                              title = gt::md("**Analysis report**"),
                              subtitle = gt::md("*Standardized residuals*"))
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
                                  title = gt::md("**Analysis report**"),
                                  subtitle = gt::md("*Moment-corrected standardized residuals*"))
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
                                  title = gt::md("**Analysis report**"),
                                  subtitle = gt::md("*Adjusted standardized residuals*"))
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


  #define the 'gt' elements for the output table of the analysis' results
  report.out <- gt::gt(report.df)

  report.out <- gt::cols_align(report.out,align="right", columns=VALUE)

  report.out <- gt::tab_header(report.out,
                           title = gt::md("**Analysis report**"),
                           subtitle = gt::md("*Chi-square and G-square test, and measures of association*"))

  report.out <- gt::cols_label(report.out,
                           STATISTIC = gt::md("**STATISTIC**"),
                           VALUE = gt::md("**VALUE**"),)

  report.out <- gt::tab_source_note(report.out,
                                           md("Thresholds for the strength of categorical association: <br>
                                           *Weak* association 0.00-0.10; *Moderate* association 0.11-0.30; *Strong* association 0.31-1.00
                                           <br> <br> after Healey J.F. 2013, *The Essentials of Statistics. A Tool for Social Research*, 3rd ed., Belmont: Wadsworth"))

  report.out <- gt::tab_options(report.out, table.font.size=tfs, source_notes.font.size=10, table.width = gt::pct(100))


  #print the various 'gt' formatted tables defined in the above steps
  print(input.table.out)
  print(exp.freq.out)
  print(chisq.values.out)
  print(relative.contrib.out)
  print(absolute.contrib.out)
  print(stand.res.out)
  print(mom.corr.stand.res.out)
  print(adj.stand.res.out)
  print(report.out)


  #create a list of results to be returned
  results <- list("crosstab"=df,
                  "exp.freq."=exp.freq,
                  "chisq.values"=chisq.values,
                  "chisq.relat.contrib."=relative.contrib,
                  "chisq.abs.contrib."=absolute.contrib,
                  "chisq.statistic"=chisq.stat,
                  "chisq.p.value"=p,
                  "chisq.p.value.perm."=p.uppertail,
                  "Gsq.statistic"=Gsquared,
                  "Gsq.p.value"=p.Gsquared,
                  "stand.resid."=stand.res,
                  "mom.corr.stand.resid."=mom.corr.stand.res,
                  "adj.stand.resid."=adj.stand.res,
                  "Phi"=phi,
                  "Phi.signed"=phi.signed,
                  "Yule's Q"=Q,
                  "Yule's Q p.value"=Q.p,
                  "Odds ratio"=or,
                  "Odds ratio CI lower boundary"=or.lower.ci,
                  "Odds ratio CI upper boundary"=or.upper.ci,
                  "Odds ratio p.value"=or.p.value,
                  "Cadj"=Cadj,
                  "Cramer's V"=V,
                  "Cramer's V CI lower boundary"=V.lower,
                  "Cramer's V CI upper boundary"=V.upper,
                  "Cramer's Vbc"=V.bc,
                  "w"=w,
                  "lambda (rows dep.)"=lambda.row.dep,
                  "lambda (cols dep.)"=lambda.col.dep,
                  "lambda (symmetric)"=lambda,
                  "tau (rows dep.)"=tau.row.dep,
                  "tau (cols dep.)"=tau.col.dep,
                  "gamma"=gamma.coeff,
                  "gamma.p.value"=gamma.p,
                  "k"=cohen.k,
                  "k CI lower boundary"=k.lower.ci,
                  "k CI upper boundary"=k.upper.ci)

  if (plot.or==TRUE) {
    visualize_odds_ratios(df, reference_level = reference_level, row_category = row_category , or.alpha = or.alpha)
  }

  return(results)
}

