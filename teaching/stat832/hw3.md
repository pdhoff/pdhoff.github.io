#### Written exercises 
Do exercises 13, 14, 16 and 18 in the course
notes on eigendecomposition (dated 1-28 or later). 


#### Heads
The dataset `heads` in the course directory
consists of six measurements on the heads of 200 soldiers in
the Swiss army. The six variables are

* MFB minimum frontal breadth
* BAM breadth of angulus mandibulae
* TFH true facial height
* LGAN length from glabella to apex nasi
* LTN length from tragion to nasion
* LTG length from tragion to gnathion

These data were obtained from the R package `Flury`.
Perform a PCA on these data.
Specifically:

1. Obtain the eigendecomposition of the sample covariance matrix, and
   interpret the first and second eigenvectors.
2. Find the correlation between the six original variables and the
   first two principal components (a six by two matrix). Make a
   scatterplot of these correlations, draw a unit circle on the plots,
   and interpret the plot.
3. For $q=0,\ldots,6$, how much variation in the data is
   explained by the best $q$-dimensional affine approximation?
4. The Swiss army is considering ordering gas masks, and are debating
   between making either two types of gas masks or four.
   If they make two types, describe how the two types should differ
   from one another. If they make four types, describe what the
   differences should be. Do you think they should make two or four types?

#### Phonemes
The `phoneme` data contains $n=4509$ periodograms, each periodogram
consisting of $p=256$ values and representing an audio sample of
someone speaking a sound. More information on these data can be found
[here](https://web.stanford.edu/~hastie/ElemStatLearn/datasets/phoneme.info.txt). 

1. Obtain the eigenvalues of the sample covariance matrix.
Make a plot of the ordered eigenvalues and the cumulative proportion
of variance explained. What fraction of the variance in the data
is explained by the first four principal components?
2. The rownames of the data matrix indicate the sounds that
each periodogram represents. Make scatterplots of the first
four principal components, indicating the different sounds by plotting
color, character or text.
3. Compute the four-dimensional vector of principal component means
for each of the five different sounds. Report these values and include
them on the previous plots.
4. Categorize each row of the data matrix to the sound category
that it is closest to, in terms of the first four principal component scores.
Make a five-by-five two-way table indicating how the different
sounds are categorized.




