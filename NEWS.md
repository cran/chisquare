# chisquare 0.8:

* Permutation-based p-value for the chi-squared statistic added;
* The calculation of the p-value for the permutation-based and Monte Carlo method has been modified from 
(1+sum(chistat.simulated>=chisq.stat)) / (B+1) to sum(chistat.simulated>=chisq.stat)/B;
* The optional plot of the permutation distribution of the chi-square statistic has been added;
* The option of plotting the permutation and the Monte Carlo distribution side-by-side has been added;
* In the plot showing the distribution of the Monte Carlo chi-squared values, the main and sub title have been amended: reference to "permutation" has been changed to "Monte Carlo";
* The minimum and the average expected count are now reported in the last rendered table;
* An annotation has been added to the first rendered table in order to give suggestions as to the selection of the most "reliable" chi-square result on the basis of different criteria (following literature; see help documentation);
* In the first rendered output table, the items featuring the annotation at the bottom of the table have been grouped and separated by an empty line to enhance readability;
* Unlike version 0.7, for determining the magnitude of the association expressed by Cadj, Cramer's V, and V biased-corrected, the coefficients are not internally converted to Cohen's w, but the Cohen's thresholds for small, medium, and large effect are dinamically internally computed on the basis of the table's df;
* For consistency, the indication of the magnitude of the association for the W coefficient is now based on the Cohen's thresholds, using the same procedure employed for Cadj, Cramer's V, and V biased-corrected;
* Scholarly literature about the average expected count added to the help documentation;
* Help documentation updated and improved content-wise.
