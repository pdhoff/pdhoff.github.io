#### Written exercises 
Do  exercises 7 and  9-12
of the normal distribution notes. 

#### Numerical exercises 
In this exercise you will examine the properties of the
sample covariance of a normal matrix under the partial
isotropy model (the "spiked covariance" model). You will
simulate data from the following model:

$$
Y \sim N_{n\times p} (0 , \Sigma\otimes I)
$$

where $\Sigma = \text{diag}( \lambda_1 + 1,\lambda_2 +1 , \lambda_3+1,1,\ldots,1)$.
For each simulation scenario, you will simulate 1000 (or more)
data matrices $Y$ compute the eigendecomposition of $Y^\top Y= V L V^\top$.
For each simulated dataset, save the values of $v^2= \text{diag}(V)^2$, the squared diagonal
elements of $V$, and the values of $L/n$.

1. Simulate data for the case $n=50$, $p=5$, and $(\lambda_1,\lambda_2,\lambda_3)=(9,6,3)$.
   Make boxplots that describe the marginal distributions of $v^2$ and the scaled
   eigenvalues.
   Repeat for the case $(\lambda_1,\lambda_2,\lambda_3)=(9,4,4)$ and compare the results.
   Explain any differences that you see.

2. Repeat 1 for the case $n=150$. Explain differences that you see.

3. Repeat 1 for the case $n=50$ and $p=50$ and interpret the results.

4. Repeat 1 for the case $n=50$ and $p=100$ and interpret the results.

