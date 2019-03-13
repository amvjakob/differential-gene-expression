process
================
Yanchao
2019-03-12

- ### Since your initial proposal, you should have decided more concretely on what methods to use for each step of your analyses, and employed some of those methods.
- ### Briefly and concisely explain your methodology and progress for the aims you have investigated so far. Which parts were modified and which parts remained the same?

#### The following are steps that have already been completed.

**Clean the data:**

We kept all the initial information from two datasets. Then we converted age to a categorical variable which has three values: larger or equal to 5-year-old, smaller than 5-year-old, and no information. The two pre-processed GEO datasets were merged after filtering the data. The combined data and merged metadata could be found in the initial analysis RMD file.

**Normalization:**

We first took log2 transformation for both datasets. Since the scales of the two datasets are different, we decided that quantile normalization was appropriate for our data after considering a few normalization methods.

**Data exploration:**

We performed some graphs (histogram and density plot) to see the distribution of the merged data. We also assigned some values for missing data and made sure there is no missing value in the combined dataset. 

**Principal Component Analysis (PCA):**

PCA was conducted on gene expression to identify the variability explained by factors available in the combined dataset. The first three principal components capture around 60% of the variability in our data. We also plotted the PC1 vs. PC2 for age, batch, and diagnosis separately. The results of these plots show that age, batch, and diagnosis are all not related to PC1 and PC2. 

**Hierarchical Clustering:**

The clustering method we applied was Average method, and our distance metric was “Euclidean”. We still want to find a way to validate our clustering result. 

**Linear regression:**

Linear regression was performed by using the “limma” package in R to identify the top differentially gene expression between control and autism cases. We identified 40 different genes, given a p-value cutoff of 0.1. Then we plotted the top 10 genes and used Hierarchical Clustering algorithms to cluster the top 40 genes that showed differential expression between control and autism cases.

#### The following steps are still in progress.



- ### What R packages or other tools are you using for your analyses? You do not need to provide your scripts in your report.

Specific packages required by our project are: cluster, pvclust, xtable, limma, GEOquery, knitr, pheatmap, stringr, ggplot2,reshape2, tidyverse.

- ### Provide the links to any markdown reports within your repo to refer to the relevant analysis.

[Initial analysis](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/Initial_analysis.md)


Rezult 

- ### What are your primary results? Were you able to answer your hypothesis? Did you have any positive results? If no, postulate a discussion as to why that may be. Provide plots and/or tables to present your results. - List some challenges that you encountered? How will you address them?

We identified 40 different genes between between control and autism cases (p-value cutoff = 0.01) by using the multiple linear regression. one of our hypothesis was that different gene expression would be detectable in comparing control and autism cases. 
