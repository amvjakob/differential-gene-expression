Technical Report
================
-------

### Introduction

Autism Spectrum Disorder (ASD) large-scale whole exome sequencing studies have found there is no single gene, but a collection of rare variants distributed across many genes that confer its manifestation[2],[3]. The variants occur in regions of multiple general transcription factors, which lend to alterations in global levels of gene expression regulation. Expression profiles have been generated with the use of RNA sequencing data from blood samples to develop predictive risk assessments of ASD in children. Previous studies have targeted peripheral blood lymphocytes (leukocyte type)[4] and whole peripheral blood (erythrocytes, thrombocytes and leukocytes)[5] to generate transcriptome signatures, yet there has been little investigation as to whether predictive models may successfully be applied to a independent datasets. We are interested in comparing the global gene expression profiles of children with ASD versus a control group from these two studies, by generating our own predictive model with machine learning methods. To further substantiate the veracity of the model, we will use the metadata to evaluate global expression signatures across age, batch and diagnosis (autism vs control). Accounting for the age variable was selected as gene regulation in individuals with ASD has been shown to fluctuate over maturation[6]. The validation of predictive risk models is fundamental to the refinement of ASD diagnostic tools and may lend to the development of gene-targeted treatments.


### Statistical Analysis

Clean the data: First, we combined data from two publicly available datasets from the Gene Expression Omnibus (GEO) database (GEO accession numbers GSE18123 and GSE25507). The second dataset was subdivided into two separate datasets, which were generated using different sequencing platforms in the study. No rubric was provided to merge these datasets by IDs successfully. For this reason, we selected the first of the two datasets (GSE25507), which uses the same ID references as our first dataset (GSE18123). We choose the interested variables such as age, batch, diagnosis from the two metadata.

We convert diagnosis to a categorical variable which has two levels(autism and control).

GSE25507 has batch information(batch 1 and batch 2); however, GSE18123 has no information about batch, so we assign "none" for the GSE18123 metadata. We want to use PCA or cluster analysis to identify whether there is a batch effect of our combined data.

There are some missing values (15 samples) for our age variable. We consider these missing points might correlate with other variables(such as gender or batch). This satisfies the condition that the missing value is missing at random. Then we perform the multiple imputations with five repetitions for our age variable. The two pre-processed GEO datasets were merged after these the step.

Normalization: We first took log2 transformation for both gene expression datasets. Since the scales of the two datasets are different, we decided that quantile normalization was appropriate for our data after considering a few normalization methods. The density plot of average gene expression value can be found in our analysis. For further details of this step, please refer to [pre\_processing\_data.Rmd]().

Afterward, we performed our exploratory analysis on the data and model the data with Limma. Several analytical approaches we applied include PCA, agglomerative hierarchical clustering, modeling with limma.

Principal Component Analysis (PCA): PCA was conducted on gene expression to identify the variability explained by factors available in the combined dataset. The first three principal components capture around 60% of the variability in our data. We also plotted the PC1 vs. PC2 for age, batch, and diagnosis separately. The results of these plots show that age, batch, and diagnosis are all not related to PC1 and PC2.

Hierarchical Clustering: The clustering method we applied was Average method, and our distance metric was “Euclidean”. We did so to determine if our variables such as age, batch, and diagnosis are well-defined clusters.

Linear regression: Linear regression was performed by using the “limma” package in R to identify the top differentially gene expression between control and autism cases. We identified 15 different genes, given a p-value cutoff of 0.1. Then we plotted the top 10 genes and used Hierarchical Clustering algorithms to cluster the top 15 genes that showed differential expression between control and autism cases.

Normality check: We need to check whether the residuals are normally distributed when we use linear regression. We use the Anderson-Darling and the Shapiro-Wilks test to check the residuals. We also us a Bonferroni correction on the significance threshold to control the number of false positives.

More details for our analysis can be found in [exploratory\_and\_limma\_analyses.md]().


### Gene Enrichment Analysis

A geneset enrichment analysis was employed to evaluate the expression values of the genes for GO term multifunctionality. Input data for this analysis was generated in the [pre-processing]( ) step, where the expression data is combined from both datasets (GEO accession numbers GSE18123 and GSE25507 (P1)), cleaned and normalised. To further prepare the data, a differential expression analysis was run by creating a design matrix with the metadata, lmfit() to fit a linear model, and ebayes() to compute logFC values and the probes IDs were annotated with their corresponding gene symbols. The input data was subsequently subsetted to gene symbol and logFC. 

Gene multifunctionality scores were inputted from the STAT540 Github repository and used to check for the degree of multifunctional bias. The method employed was a Spearman’s correlation, which outputted a coefficient value close to zero, indicating little to no multifunctional bias. 

Finally, a Precision-Recall method was utilised to calculate the p-value of GO terms after a multifunctionality correction. This enrichment method put emphasis on highly ranked genes, identifying enrichment that could potentially impact the weight of the predictor genes generated from the machine learning analysis. GO.xml file was accessed using ErmineR to place genes into sets.  As expected from the Spearman’s correlation coefficient, GO terms were only minimally corrected in the enrichment analysis, with the largest log p-value adjustment at 0.4714262. 

These results overall indicate a low level of multifunctional bias. Furthermore, the top ranked genesets for bias did not contain any of the predictor genes generated in the machine learning analysis. 
