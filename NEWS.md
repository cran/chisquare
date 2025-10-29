# What's new in chisquare version 1.2

## Statistical methodology improvements

* Exact p-value calculation for permutation and Monte Carlo tests: Implemented the Phipson & Smyth (2010) method for computing exact p-values when permutations are randomly drawn. The formula p = (B + 1)/(m + 1) ensures that p-values are never zero and that Type I error rates are correctly controlled. This supersedes the previous implementation that could underestimate p-values and inflate Type I error rates;

## New association measures and effect sizes

* Maximum-corrected version of the C coefficient added;
* 95 percent confidence interval around Goodman-Kruskal lambda (asymmetric and symmetric) added;
* 95 percent confidence interval around Goodman-Kruskal tau added;
* Enhanced effect size interpretation for chi-square-based measures: uncorrected measures (Phi, C, Cramer's V) are now assessed against Cohen's thresholds adjusted by their maximum achievable values, while corrected measures use standard Cohen's thresholds;
* Added a detailed effect size interpretation thresholds table showing exact threshold values for all association measures;

## New residual analysis methods

* Goodman-Kruskal residuals added;
* PEM (Percentage of Maximum Deviation from Independence) added;
* Standardised median polish residuals and adjusted standardised median polish residuals added, providing robust alternatives for cell residual analysis that are resistant to masking and swamping effects from outliers;

## User interface and output improvements

* Added a new `gt.tables` component to the returned results that contains all formatted gt tables; this allows the tables to be rendered any time without needing to re-run the analysis;
* Modification to the code to avoid automatic rendering of all tables in the R console, improving execution speed and preventing browser overload; tables are now returned in a structured list (see previous bullet point) but only rendered when setting the new `render.all.tables` parameter to TRUE;

## Documentation improvements

* Major restructuring of the help documentation for improved clarity and pedagogical flow: content is now organised into six thematic sections (Hypothesis Testing and Power Analysis; Decomposition of the Chi-Square Statistic; Table Standardisation; Chi-Square-Maximising Table; Post-Hoc Analysis; Measures of Categorical Association);
* Table standardisation methodology (Iterative Proportional Fitting) now documented separately from specific coefficient applications, enhancing conceptual clarity;
* Association measures reorganised by statistical foundation: Chi-square-based measures, Margin-free measures, and PRE (Proportional Reduction in Error) measures, facilitating coefficient selection based on analytical needs;

## Bug fixes and minor improvements

* Annotation at the bottom of the last rendered table removed since no longer relevant;
* Typo in the title of the absolute contribution to the chi-square statistic table fixed;
* Updates to the help documentation's content and layout;
* Bibliographical references updated.
