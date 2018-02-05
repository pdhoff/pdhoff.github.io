#### Written exercises 
Do  exercises 1 - 6  from the SVD notes. 

#### AMMI
Do an AMMI analysis of the `exports` dataset in the course directory. 
These data give yearly export totals reported to UN Comtrade for 
131 countries from 2001 through 2010, given in 2010 inflation 
adjusted dollars. Transform these data to the log scale before 
doing any analysis. 

1. Do a two-way ANOVA, give the ANOVA table, and report the outcome 
of the two $F$-tests for 
non-zero additive row effects and non-zero additive column effects. 

2. Obtain the SVD of the data matrix after subtracting off the 
row and column means. Make a scree plot of the singular values, and make a
table of the proportion of variance explained by the best low-rank 
matrix approximations of rank zero to nine, and comment. Why don't 
we need to consider a rank-10 approximation?

3. Fit a rank-2 AMMI model, and make a plot of the estimated multiplicative 
effects, using meaningful plotting labels. Interpret the plot. 




