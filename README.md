# chisquare 0.8

vers 0.7:
* Indication of effect size added to Phi, Phi signed, Phi corrected, Cadjusted;
* Changes and fix in the method underlying the indication of the effect size for Cramer's V and bias-corrected V;
these changes (following Sheskin 2011) allow to provide indication of effect size for tables of any size;
* Wald 95 percent confidence interval around the Monte Carlo (permutation-based) chi-square p-value added;
* The 95 percent CI around the permuted p-value is reported in the annotations at the bottom of the first produced table, and in the `Analysis report` table;
* Phi-corrected added;
* W coefficient added (including magnitude of effect size and bootstrap 95% CI);
* W is also reported in the annotation at the bottom of the first tables produced by the function;
* Corrected version of Goodman-Kruskal's lambda (both aymmetric and asymmetric) added;
* In the annotations at the bottom of the first table produced by the function, Goodman-Kruskal symmetric lambda has been replaced by its corrected version for its ability to addresses skewed or unbalanced marginal probabilities which create problems to the traditional lambda;
* The magnitude of the effect indicated by the odds ratio (introduced in ver. 0.6) is not returned anymore, but only reported in the `Analysis report` table for consistency with other effect size magnitude indication;
* Labels used for the magnitude of effect size slightly modified to enhance clarity and consistency;
* The returned results are organised into a nested list format for better clarity and accessibility: `input.table`, `chi.sq.related.results`, `G.square`, `residuals`, `chi.sq.based.assoc.measures`, and `non.chi.sq.based.assoc.measures`;
* Help documentation updated and improved content-wise;
* An error in the help documentation has been amended: the adjusted chi-square statistic was described as produced by "dividing" the chi-square value by (N-1)/N; "dividing" was meant to read "multiplying";
* Relevant scholarly references added.


vers 0.6
* chi-square test adjusted by the (N-1)/N correction added;
* magnitude of the effect size as indicated by the odds ratio added;
* help documentation updated, and improved content-wise;
* relevant scholarly references added.

vers 0.5
* the plot of odds ratios is now returned for cross-tabs featuring 2 rows and at least 2 columns;
* to enhance consistency and clarity, the parameter 'row_category' has been renamed to 'row.level'; the parameter 'reference_level' has been renamed to 'reference.level'; 
* Haldane-Anscombe correction added in the calculation of odds ratio;
* to improve clarity, minor adjustments have been carried out to the error messages returned when the parameter 'plot.or' is 'TRUE' but the input cross-tab does not meet the needed requirements (2xk, where k >= 2);
* help documentation updated and improved content- and layout-wise.

vers 0.4
* the facility to plot the odds ratios and 95% confidence intervals for table of size 2xk (where k>2) has been added;
* the above added facility rests on a new internal function (visualize_odds_ratios);
* amendments and updates to the help documentation.

vers 0.3
* Moment-corrected standardized residuals added;
* 95% confidence interval for Cramer's V added; 
* The output version of the input cross-tabulation now reports some essential analytical results in an annotation at the bottom of the table;
* Fix to the output version of the input cross-tabulation when cells feature a frequency of 0;
* The Goodman-Kruskal's tau reported in previous versions is now indicated as "columns dependent"", and the "row dependent" version of the same coefficent is now calculated and reported.

vers 0.2
* Relative contributions to the chi-square added;
* Contributions to the chi-square previously reported in ver_0.1 are now called absolute contribution;
* Bibliographical reference added for the abovementioned contributions;
* Significance level of the reported standardized and adjusted standardized residuals can be corrected using
the Sidak's alpha adjustment method;
* G-square test added;
* Symmetric Goodman-Kruskal's lambda and Goodman-Kruskal's tau added;
* Cohen's k (with 95% CI) added;
* Odds ratio (with p value and 95% CI) added (for 2x2 tables);
* Effect size for Cramer's V and Cramer's V bias-corrected added to the table of the test's results; returned for tables with up to 5 degrees of freedom according to Cohen 1988;
* Degrees of freedom added to the output table reporting the test's result;
* The distribution of the permuted chi-square statistic can be optionally plotted;
* Fix in the calculation of the cells contribution (in percent) to the chi-square;
* Output tables formatted using the 'gt' package instead of 'formattable';
* Option added to customize the font size of the output tables;
* Improvements to the content and layout of the help documentation.


vers 0.1
first release to CRAN
