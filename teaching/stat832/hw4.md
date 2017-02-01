#### Written exercises 
Do  exercise 1 from the LDA notes and exercises 1 and 2 from the normal 
distribution notes. 

#### Training and test error in LDA 
Compare some classifiers for the `phoneme` data, and evaluate 
their out of sample predictive performance. Randomly divide the $n=4509$ 
periodograms into two groups of (nearly) equal size, one to train your 
classifiers and one to test them. Your first type of classifier is the 
one from the Homework 2, and is obtained as follows:

1. Using the training data, obtain the matrix $V$ that transforms the 
   training data to the principal components. 
2. For case $i$ in the training set 
   let $x_i$ be the vector of the first $r$ principal components. 
   Compute sample means of the $x_i$'s for each class. Compute the 
   in-sample missclassification rate based on classifying each $i$ to 
   the class whose group mean is closest to $x_i$. 
3. Rotate the test data by multiplying by the matrix $V$ obtained from 
   the training data in part 1. Compute the out-of-sample 
   missclassification rate based on classifying each test case $i$ to 
   the class whose group mean (computed in 2 from the training data) 
   is closest to $x_i$.  
4. Do 1-3  for each $r$ from 1 to some number where it 
   seems pointless to continue. At the end, you wil have two 
   misclassification curves, one being in-sample and the other 
   out-of-sample. 

Now repeat the above procedure but using LDA to obtain your classification 
rule. Instead of using the first $r$ PCs, use the first $r$ linear 
discriminants. Compare the missclassification rates here to using 
the PCA approach, and comment.  



