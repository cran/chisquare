#' R function for chi-square test of independence
#'
#' The function performs the chi-square test of independence on the input contingency table,
#' calculates various measures of categorical association, and returns standardized and adjusted
#' standardized residuals. The p value associated to the chi-square statistic is also calculated on the basis
#' of a permutation-based procedure. Nicely-formatted output tables are rendered.
#'
#' The following \strong{measures of categorical associations} are produced by the function:
#'  \itemize{
#'   \item Phi (only for 2x2 tables)
#'   \item Phi signed (only for 2x2 tables)
#'   \item Yule's Q (only for 2x2 tables)
#'   \item Adjusted contingency coefficient C
#'   \item Cramer's V
#'   \item Bias-corrected Cramer's V
#'   \item Cohen's w
#'   \item Goodman-Kruskal's lambda (asymmetric)
#'   \item Goodman-Kruskal's gamma
#' }
#'
#' The \strong{p value} of the observed chi-square statistic is also calculated on the basis of a permutation-based approach,
#' using B random tables created under the Null Hypothesis of independence. For the rationale of this approach,
#' see for instance the description in Beh E.J., Lombardo R. 2014, Correspondence Analysis: Theory, Practice
#' and New Strategies, Chichester, Wiley, pages 62-64.\cr
#'
#' The \strong{randomized p value} is calculated as follows: \cr
#' \eqn{(1 + sum (chistat.perm > chisq.stat)) / (1 + B)}, where
#' \emph{chistat.perm} is a vector storing the B chi-square statistics generated under the Null Hypothesis, and
#' \emph{chisq.stat} is the observed chi-square statistic. For the logic of the calculation, see for example
#' Baddeley et al., "Spatial Point Patterns. Methodology and Applications with R", CRC Press 2016: 387.\cr
#'
#' The \strong{adjusted standardized residuals} are calculated as follows: \cr
#' \eqn{stand.res[i,j] / sqrt((1-sr[i]/n)*(1-sc[j]/n))}, where \emph{stand.res} is the standardized residual for cell \emph{ij},
#' \emph{sr} is the row sum for row \emph{i}, \emph{sc} is the column sum for column \emph{j}, and
#' \emph{n} is the table grand total. The \emph{adjusted standardized residuals} may prove useful since it has
#' been demonstrated that the standardized residuals tend to underestimate the significance of differences in small samples.
#' The adjusted standardized residuals correct that deficiency.\cr
#'
#' The \strong{bias-corrected Cramer's V} is based on Bergsma, Wicher (2013).
#' "A bias correction for Cramér's V and Tschuprow's T". Journal of the Korean Statistical Society. 42 (3): 323–328,
#' https://doi.org/10.1016%2Fj.jkss.2012.10.002\cr
#'
#' For the \strong{other measures of categorical association} provided by the function, see for example
#' Sheskin, D. J. (2011). Handbook of Parametric and Nonparametric Statistical Procedures, Fifth Edition (5th ed.). Chapman and Hall/CRC: 675-679,
#' 1415-1427.\cr
#'
#' \strong{Note} that the \emph{Phi} coefficient is based on the chi-square statistic as per Sheskin's equation 16.21, whereas the
#' \emph{Phi signed} is after Sheskin's equation 16.20. The \emph{2-tailed} p value of \emph{Goodman-Kruskal's gamma} is based on the
#' associated z-score calculated as per Sheskin's equation 32.2.
#'
#'
#' @param data Dataframe containing the input contingency table.
#' @param B  Number of simulated tables to be used to calculate the Monte Carlo-based p value (999 by default).
#' @param format Takes \emph{short} (default) if the dataset is a dataframe storing a contingency table; if the
#' input dataset is a dataframe storing two columns that list the levels of the two categorical variables,
#' \emph{long} will preliminarily cross-tabulate the levels of the categorical variable in the 1st column against
#' the levels of the variable stored in the 2nd column.
#'
#' @return The function produces a number of \strong{output tables}, which are nicely formatted with the help of the
#' \code{\link[formattable]{formattable}} function out from the the \emph{formattable} package.
#' The output tables are listed below:
#'
#'  \itemize{
#'   \item Input contingency table
#'   \item Expected frequencies
#'   \item Cells' chi-square value
#'   \item Cells' contribution (in percent) to the chi-square statistic (cells in blue feature a larger-than-average
#'   contribution)
#'   \item Standardized residuals (RED for residuals smaller than -1.96, GREEN for residuals larger than +1.96)
#'   \item Adjusted standardized residuals (colour same as above)
#'   \item Table of output statistics, p values, and association measures
#' }
#'
#' Also, the function returns a \strong{list containing the following elements}:
#'  \itemize{
#'   \item \emph{crosstab}: input contingency table
#'   \item \emph{exp.freq.}: table of expected frequencies
#'   \item \emph{chisq.values}: cells' chi-square value
#'   \item \emph{chisq.contrib.}: cells' contribution (in percent) to the chi-square statistic
#'   \item \emph{chisq.statistic}: observed chi-square value
#'   \item \emph{chisq.p.value}: p value
#'   \item \emph{chisq.p.value.perm.}: p value based on B permuted tables
#'   \item \emph{stand.resid.}: table of standardized residuals
#'   \item \emph{adj.stand.resid.}: table of adjusted standardized residuals
#'   \item \emph{Phi}: Phi coefficient (only for 2x2 tables)
#'   \item \emph{Phi signed}: signed Phi coefficient (only for 2x2 tables)
#'   \item \emph{Yule's Q}: Q coefficient (only for 2x2 tables)
#'   \item \emph{Cadj}: adjusted contigency coefficient C
#'   \item \emph{Cramer's V}: Cramer's V coefficient
#'   \item \emph{Cramer's Vbc}: bias-corrected Cramer's V coefficient
#'   \item \emph{w}: Cohen's w
#'   \item \emph{lambda (rows dep.)}: Goodman-Kruskal's lambda coefficient (considering the rows being the dependent variable)
#'   \item \emph{lambda (cols dep.)}: Goodman-Kruskal's lambda coefficient (considering the columns being the dependent variable)
#'   \item \emph{gamma}: Goodman-Kruskal's gamma coefficient
#'   \item \emph{gamma.p.value}: 2-sided p value for the Goodman-Kruskal's gamma coefficient
#' }
#'
#' The following examples, using in-built datasets, can be run to fimiliarise with the function:\cr
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
#' -perform the test, specifying that the input table is in 'long' format\cr
#' result <- chisquare(mytable, format="long")\cr
#'
#'
#' @keywords chiperm
#'
#' @export
#'
#' @importFrom formattable formatter formattable style
#' @importFrom stats pchisq addmargins r2dtable pnorm
#'
#' @examples
#' #perform the test on the in-built 'diseases' dataset
#' result <- chisquare(diseases, B=99)
#'
#'
chisquare <- function(data, B=999, format="short"){
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
  chisq.values <- round(calc(df)$chisq.values,2)
  exp.freq <- round(calc(df)$exp.freq,2)

  #p value
  p <- round(pchisq(q=chisq.stat, df=(nr-1)*(nc-1), lower.tail=FALSE), 3)

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

  #cells' contribution to the chi-sq
  contrib.to.chisq <- round(chisq.values/chisq.stat * 100,2)

  #average cell contribution to chi-sq statistic
  #to be used as threshold to give different colours
  #in the table produced later on with formattable
  average.contr <- 100 / (nr*nc)

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
  } else {
    Q <- "-"
  }

  #adjusted contingency coeffcient
  C <- sqrt(chisq.stat / (n + chisq.stat))
  Cmax <- sqrt((min(nr,nc)-1) / min(nr,nc))
  Cadj <- round(C / Cmax,3)

  #Cramer's V
  V <- round(sqrt(chisq.stat / (n * min(nr-1, nc-1))), 3)

  #bias-corrected Cramer's V
  phi.new <- max(0, (chisq.stat / n) - ((nc-1)*(nr-1)/(n-1)))
  k.new  <-  nc-((nc-1)^2 / (n-1))
  r.new  <-  nr-((nr-1)^2 / (n-1))
  V.bc <- round(sqrt(phi.new / min(k.new-1, r.new-1)), 3)

  #Cohen's w
  w <- round(V * sqrt(min(nr, nc)-1),3)
  w_label <- ifelse(w <= 0.3, "small effect",
                    ifelse(w <= 0.5, "medium effect",
                           "large effect"))

  #Goodman-Kruskal's lambda
  E1 <- n - max(sr)
  E2 <- sum(sc - apply(df, 2, max))
  lambda.row.dep <- round((E1-E2)/E1,3)

  E1 <- n - max(sc)
  E2 <- sum(sr - apply(df, 1, max))
  lambda.col.dep <- round((E1-E2)/E1,3)

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
  gamma.p <- calc.gamma(df)$GK.gamma.p

  gamma.p.to.report <- ifelse(gamma.p < 0.001, "< 0.001",
                        ifelse(gamma.p < 0.01, "< 0.01",
                               ifelse(gamma.p < 0.05, "< 0.05",
                                      round(gamma.p, 3))))

  #standardized residuals
  stand.res <- round((df - exp.freq) / sqrt(exp.freq), 2)

  #create a copy of the previous dataframe to be populated
  #with the values of the adjusted standardized residuals
  adj.stand.res <- stand.res

  #adjusted standardized residuals
  for (i in 1:nr) {
    for (j in 1:nc) {
      adj.stand.res[i,j] <- round(stand.res[i,j] / sqrt((1-sr[i]/n)*(1-sc[j]/n)), 2)
    }
  }

  #define a function to conditionally format the cells' contribution to the chi-sq;
  #the function will be used in a subsequent'formattable' function
  myForm1 <- function(x) {
    formattable::formatter('span', style = x ~ style(color = ifelse(x > average.contr, "blue", "black")))
  }

  #define a function to conditionally format the standardized residual values;
  #the function will be used in a subsequent'formattable' function
  myForm2 <- function(x) {
    formattable::formatter('span', style = x ~ style(color = ifelse(x > 1.96, "green", ifelse(x < -1.96, "red", "black"))))
  }


  #define the outputs for the formatted tables
  df.out <- formattable::formattable(as.data.frame(addmargins(as.matrix(df))), align = "c")

  exp.freq.out <- formattable::formattable(as.data.frame(exp.freq), align = "c")

  chisq.values.output <- formattable::formattable(round(chisq.values, 3), align = "c")

  contrib.to.chisq.output <- formattable::formattable(contrib.to.chisq, lapply(contrib.to.chisq, myForm1), align = "c")

  stand.res.output <- formattable::formattable(stand.res, lapply(adj.stand.res, myForm2), align = "c")

  adj.stand.res.output <- formattable::formattable(adj.stand.res, lapply(adj.stand.res, myForm2), align = "c")

  #define a vector for the labels of the statistics to be returned
  statistics <- c("Chi-square",
                  "Phi",
                  "Phi signed",
                  "Yule's Q",
                  "Cadj",
                  "Cramer's V",
                  "Cramer's V bias-corrected",
                  paste0("Cohen's w", " (", w_label, ")"),
                  "Goodman-Kruskal's lambda (rows dependent)",
                  "Goodman-Kruskal's lambda (columns dependent)",
                  "Goodmak-Kruskal's gamma")

  #define a vector for the statistics to be returned
  values.to.report <-c(paste0(chisq.stat, " (p: ", p.to.report, "; permuted p: ", p.uppertail.to.report, ")"),
                       phi,
                       phi.signed,
                       Q,
                       Cadj,
                       V,
                       V.bc,
                       w,
                       lambda.row.dep,
                       lambda.col.dep,
                       paste0(gamma.coeff, " (p: ", gamma.p.to.report, ")"))

  #create a dataframe storing the two above-defined vectors
  report.df <- data.frame(STATISTIC=statistics, VALUE=values.to.report)

  #create an object that stores the formatted dataframe created at the previous step
  report.df.out <- formattable::formattable(report.df, align = c("l", "r"))

  #print the various formatted tables
  print(df.out)
  print(exp.freq.out)
  print(chisq.values.output)
  print(contrib.to.chisq.output)
  print(stand.res.output)
  print(adj.stand.res.output)
  print(report.df.out)

  results <- list("crosstab"=df,
                  "exp.freq."=exp.freq,
                  "chisq.values"=chisq.values,
                  "chisq.contrib."=contrib.to.chisq,
                  "chisq.statistic"=chisq.stat,
                  "chisq.p.value"=p,
                  "chisq.p.value.perm."=p.uppertail,
                  "stand.resid."=stand.res,
                  "adj.stand.resid."=adj.stand.res,
                  "Phi"=phi,
                  "Phi.signed"=phi.signed,
                  "Yule's Q"=Q,
                  "Cadj"=Cadj,
                  "Cramer's V"=V,
                  "Cramer's Vbc"=V.bc,
                  "w"=w,
                  "lambda (rows dep.)"=lambda.row.dep,
                  "lambda (cols dep.)"=lambda.col.dep,
                  "gamma"=gamma.coeff,
                  "gamma.p.value"=gamma.p)

  return(results)
}
