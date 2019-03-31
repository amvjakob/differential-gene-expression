Untitled
================
Yanchao
2019-03-28

Clean and organize the data
---------------------------

Clean the data: First, we combined data from two publicly available datasets from the Gene Expression Omnibus (GEO) database (GEO accession numbers GSE18123 and GSE25507). The second dataset was subdivided into two separate datasets, which were generated using different sequencing platforms in the study. No rubric was provided to merge these datasets by IDs successfully. For this reason, we selected the first of the two datasets (GSE25507), which uses the same ID references as our first dataset (GSE18123). We choose the interested variables such as age, batch, diagnosis from the two metadata.

We convert diagnosis to a categorical variable which has two levels(autism and control).

GSE25507 has batch information(batch 1 and batch 2); however, GSE18123 has no information about batch, so we assign "none" for the GSE18123 metadata. We want to use PCA or cluster analysis to identify whether there is a batch effect of our combined data.

There are some missing values (15 samples) for our age variable. We consider these missing points might correlate with other variables(such as gender or batch). This satisfies the condition that the missing value is missing at random. Then we perform the multiple imputations with five repetitions for our age variable. The two pre-processed GEO datasets were merged after these the step.

Normalization: We first took log2 transformation for both gene expression datasets. Since the scales of the two datasets are different, we decided that quantile normalization was appropriate for our data after considering a few normalization methods. The density plot of average gene expression value can be found in our analysis. For further details of this step, please refer to [pre\_processing\_data.Rmd](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/Pre_processing.md).

Statistical methods
-------------------

Afterward, we performed our exploratory analysis on the data and model the data with Limma. Several analytical approaches we applied include PCA, agglomerative hierarchical clustering, modeling with limma.

Principal Component Analysis (PCA): PCA was conducted on gene expression to identify the variability explained by factors available in the combined dataset. The first three principal components capture around 60% of the variability in our data. We also plotted the PC1 vs. PC2 for age, batch, and diagnosis separately. The results of these plots show that age, batch, and diagnosis are all not related to PC1 and PC2.

Hierarchical Clustering: The clustering method we applied was Average method, and our distance metric was “Euclidean”. We did so to determine if our variables such as age, batch, and diagnosis are well-defined clusters.

Linear regression: Linear regression was performed by using the “limma” package in R to identify the top differentially gene expression between control and autism cases. We identified 15 different genes, given a p-value cutoff of 0.1. Then we plotted the top 10 genes and used Hierarchical Clustering algorithms to cluster the top 15 genes that showed differential expression between control and autism cases.

Normality check: We need to check whether the residuals are normally distributed when we use linear regression. We use the Anderson-Darling and the Shapiro-Wilks test to check the residuals. We also us a Bonferroni correction on the significance threshold to control the number of false positives.

More details for our analysis can be found in [exploratory\_and\_limma\_analyses.md](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/exploratory_and_Limma_analysis.md).

Statistical result
------------------

We used Limma to identify the top differentially genes To do so we fit each gene with a linear regression model as follows:
**Y**<sub>gene</sub> = *β*<sub>0</sub> + *β*<sub>1</sub>**a****g****e** + *β*<sub>2</sub>**b****a****t****c****h** + *β*<sub>3</sub>**d****i****a****g****n****o****s****i****s** + *ϵ*

where,

**Y**<sub>gene</sub> represents the vector of expression levels for each gene.

*β*<sub>0</sub> is the intercept vector.

**a****g****e** is a continuous variable that indicates people age between 1 and 17.5. *β*<sub>1</sub> is the coefficient of the **c****e****l****l****\_****t****y****p****e**.

**b****a****t****c****h** is a binary variable with 3 levels(batch 1, batch 2 and none). *β*<sub>2</sub> is the coefficient of the **b****a****t****c****h**.

**d****i****a****g****n****o****s****i****s** is a binary variable with two levels (autism and control). *β*<sub>3</sub> is the coefficient of the **d****i****a****g****n****o****s****i****s**.

*ϵ* is a variable of residuals for a certain gene.

We identified 15 different genes between control and autism cases (p-value cutoff = 0.01) by using the multiple linear regression. One of our hypotheses was that different gene expression would be detectable in comparing control and autism cases, concerning age and batch. Additionally, our PCA results show that age, diagnosis (control and autism) and batch are all not related to first two of the variability we are observing in our data.
