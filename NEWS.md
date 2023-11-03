# chisquare 0.7:

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
