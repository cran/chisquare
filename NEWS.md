# What's new in chisquare version 1.0:

* Chi-square-maximising table added; it is returned and rendered;
* Ancillary function to compute the Chi-square-maximising table added;
* Chi-squared value computed on the Chi-square-maximising table is returned and reported in the annotation of the
rendered chi-square-maximising table; it is also reported in the last rendered table;
* Cramer V corrected added; it is returned and reported in the first and last rendered tables;
* Phi.max in now returned and reported in the last rendered table;
* Changes and fixes to the ancillary function for the computation of phi.max, which was producing a wrong value in some instances;
* Standardised version of the input table added; it is returned and rendered;
* Internal ancillary function for the computation of the standardised tables added (based on Iterative Proportional Fitting);
* New parameters added to allow the user to select the type of target marginals to be used for table standardisation;
* Cramer's V computed on the standardised table added;
* 1-(Cramer's V/Cramer's V standardised) added;
* Verbal articulation of the strength of the association measured by Yule's Q added;
* Yule's Y, its p-value, and the verbal articulation of the strength of the association added;
* Independent odds ratio returned and rendered for tables larger than 2x2;
* Odds ratio, Yule's Q, and Yule's Y added to annotations in the first rendered table;
* Quetelet index and IJ association factor added: they are returned and also rendered in two dedicated tables;
* Table of adjusted standardised counts added: it is returned and rendered in a dedicated table;
* The object "residuals" contained in the returned list has been renamed as "post.hoc";
* In the first rendered table, the report of the symmetric corrected Goodman-Kruskal's lambda replaced by the report of the (uncorrected) lambdas;
* Minor modification in the style of the annotations in the first rendered table: Italics has been dropped to enhance readability. It has been kept in the last part of the annotation, where suggestions as to the "best" chi-squared analysis are given;
* Help documentation updated content-wise, with some layout restructing.
