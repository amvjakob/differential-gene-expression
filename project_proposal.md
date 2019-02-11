# Project proposal - Y0ung-parents

## Motivation and background work

The most common neurodevelopmental disorders, Intellectual Disability (ID) and Autsim Spectrum Disorder (ASD), arise from three causative factors: genes, epigenetics and the environment [1]. Large-scale whole exome sequencing studies have found there is no single gene, but a collection of rare variants distributed across many genes that confer the manifestation of ASD [2, 3]. The variants occur in regions of multiple general transcription factors, which lend to alterations in global levels of gene expression regulation. Expression profiles have been determined with the use of RNA sequencing data from blood samples to develop predictive risk assessments of ASD in children.  Previous studies have targeted peripheral blood lymphocytes [4] and whole blood [5] to generate transcriptome signatures, yet there has been no investigation as to whether these models may be cross-validated.  We are interested in comparing the global gene expression profiles of children with ASD from these two studies, while also utilising machine learning to evaluate the robustness of each data set as true predictive models. To further substantiate the veracity of the models, we will use the metadata to evaluate global expression signatures across age both intra- and inter-data sets.  Accounting for this variable is necessary as gene regulation in individuals with ASD has been shown to fluctuate over maturation [6]. The validation of predictive risk models is fundamental to the refinement of ASD diagnostic tools and may lend to the development of gene-targeted treatments.

[1](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4185273/?_escaped_fragment_=po=0.241546), [2](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4185273/?_escaped_fragment_=po=0.241546#R5),  [3](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4185273/?_escaped_fragment_=po=0.241546#R116),  [4](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3040743/), [5](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0049475), [6](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1002592)


## Division of labour 

| Name | Background | Affiliation | Responsibilites | 
| ------------- | ------------- | ------------- | ------------- |
| Benson Chang | Biotechnology, Biological Engineering | Genome Science and Technology |  Analysis literature validation, Data visualization, Statistical analysis, Data interpretation |
| Virginia Pichler |  Genomic Epidemiology| Microbiology and Immunology | Literature search, Project planning, Data inspection/Quality control |
| Yanchao Luo | Statistical methods and Empirical Likelihood | Department of Statistics | Data cleaning & processing, Statistical analysis (hypothesis testing), Data visualization | 
| Anthony Jakob | Biological Engineering | Biomedical Engineering | Machine learning model |

## Dataset

### General description and characteristics of the data

We are working with two different datasets, both data from microarrays for expression profiling.
In both datasets, RNA was amplified, labeled and hybridized, the source of the RNA being peripheral blood lymphocytes for the first dataset and peripheral blood cells for the second dataset.
In both datasets, the subjects are children suffering from Autism Spectrum Disorder (ASD) or not.

In both datasets, the samples assign a measured expression value to the corresponding gene reference ID.
In the first sample, the measured value is the RMA signal intensity.
In the second sample, the measured value is the Probe Logarithmic Intensity Error (PLIER) signal intensity.
The samples from both datasets each contain 54613 rows and the same gene reference ids, as the same kits were used.


### Technology used to generate the data.
#### Dataset 1 (paternal): 
- Expression profiling by array
- Qiagen Qiaquick kit used on blood draw
- Double round amplification, followed by biotin-labelling using Affymetrix's GeneCHip Two-Cyle Target Labeling kit
- Checks for evenly distributed range of transcript size, verification of fragmentation
- Hybridization cocktails hybridized on Affymetrix Human Genome U133 Plus 2.0 Array (in situ oligonucleotide)
 - Washing
 - Scanned on Affymetrix's GeneChip Scanner 3000 7G
- Extraction of raw signal intensity from scanned images of the array.
- The scanning was performed on two different machines, however 
> Gene expression levels were not adjusted for possible batch effects as algorithms that attempt to adjust for batch effects also alter the gene expression distribution.
- Analysis of covariance of batch numbers.
- MAS 5.0 was used for the analysis of gene expression distribution (does not alter it)
- RMA (uses quantile normalization and my remove group level differences in gene expression distribution) was used for the gene expression analysis "looking for specific gene expression differences between groups."


#### Dataset 2 (predictive):
- Expression profiling by array
- Trizol extraction of total RNA according to manufacturer
- Generation of biotin-labeled cRNA according to Affymetrix protocols
- Quantification (A260) and fragmentation of cRNA
- Hybridization of fragmented cRNA on GeneChips (Affymetrix Human Genome U133 Plus 2.0 Array)
- Scanning using Affymetrix GeneChip scanner 3000 at 2.5 microm resultion
- Affymetrix Human Gene 1.0 ST Array was also used
- The final recored signal intensity is the PLIER of the samples

### Description of data
#### Dataset 1:
- 663.4 MB
- 146 samples
- 54613 rows (genes) per sample
- 2 columns (ID_REF, VALUE) per sample
- Age of child and parents
- (No gender)
- Scan batch number

#### Dataset 2:
- 1.2 GB
- 285 samples
- 54612 rows(genes) per sample
- 2 columns (ID_REF, VALUE) per sample
- Gender
- Age of child
- Ethnicity of child
- Diagnosis
- Ethnicity
- Other diseases

## Aims and methodology

We are interested in the differential gene expression of patients with/without autism. We are also interested in how the gene expression profiles of patients with/without autism differ at certain age groups. we have found two datasets [1](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25507), [2](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18123) from two different studies that explore highly similar variables (gene expression of 54613 genes, age, control vs. autism).

The [1st](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE25507) dataset is limited in information compared to the [2nd](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE18123) dataset, so the bulk of our primary analysis will be performed on the 2nd dataset, while the 1st dataset can be used to verify the predictive modelling generated from the 2nd dataset. We will use normalization methods to ensure our comparative analysis of both datasets is consistent. The primary study from the 2nd dataset has suggested that certain genes are significantly correlated with disorders associated with autism, we can analyze and verify whether this set of genes is also observed in the sample patients with autism from the first dataset. We will conduct PCA on gene expression profiles in the 2nd dataset to identify variability as well as to explore possible batch effect of the data since 2nd dataset contains information on sample batches. We can apply clustering in the gene expression samples to see if the identified groups correspond to control patients and patients with autism.

Both datasets did not look at the significance of sample patient age affecting gene expression profiles, even though the 2nd dataset acknowledged that patient age at blood draw significantly influenced different gene expression levels of certain genes. We could rebuild a prediction model based on the 2nd dataset, but also factor in subgrouping samples into different age groups. Statistical methods involved in building our new prediction model can include linear regression and logistic regression (because our response will only have two values - control vs autism). Building our new prediction model can also involve machine learning tools (ie. random forest algorithm) to help us verify the results of our prior clustering analysis. Our prediction model can also be verified by predicting control vs autism of the samples in the 1st dataset to test the sensitivity and accuracy of our model.

@vpichler
