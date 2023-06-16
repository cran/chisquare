# chisquare 0.5:

* the plot of odds ratios is now returned for cross-tabs featuring 2 rows and at least 2 columns;
* to enhance consistency and clarity, the parameter 'row_category' has been renamed to 'row.level'; the parameter 'reference_level' has been renamed to 'reference.level'; 
* Haldane-Anscombe correction added in the calculation of odds ratio;
* to improve clarity, minor adjustments have been carried out to the error messages returned when the parameter 'plot.or' is 'TRUE' but the input cross-tab does not meet the needed requirements (2xk, where k >= 2);
* help documentation updated and improved content- and layout-wise.
