# ML models results

All models were fitted using Python and the library scikit-learn.

## Datasets
All analyses were performed on every dataset individually. Training was done on 70% of the samples. Since the second dataset itself contains a training and a test set, those were combined before 70% was used for training.

## Procedure
The expression data was extracted from the series matrix txt files.
Lowly expressed genes (expression value under 150 in at least one sample) were eliminated from the dataset.
The data was centered and normalized, and the rows shuffled.
Different models were then trained using the data, the idea being to use a L1-penalty to let the model select "relevant" features, and to compare the selected genes across the different models, datasets and types of analyses (compare with results from statistic analysis).

### SGDClassifier (Stochastic Gradient Descent)
- Cross-validation to find optimal C.
- Hinge loss, elastic net penalty (mixture of L1 and L2), L1-ratio 0.3, maximum of 100 iterations

### LogisticRegression
- Cross-validation to find optimal C.
- L1-penalty

### RandomForestClassifier
- 100 random trees with a maximal depth of 20

## Results (so far)
82.1% classification accuracy using 82 genes (SGDClassifier).

96.6% classification accuracy using 40 genes (LogisticRegression). This result needs investigation.

90.0% classification accuracy using 3949 genes (RandomForestClassifier). We need to train this model using already selected features.

## Next steps
- Train model on merged datasets
- Implement a pipeline that first selects "relevant" genes without removing lowly expressed genes and then fits a second model on top of it (e.g. random forest classifier that does not perform feature selection on its own).
- Compare the selected genes across models and datasets (e.g. use one dataset to train and the other dataset to validate the model's accuracy).
- Compare the genes selected by the models to the statistically relevant genes.
- Add age and batch number as features
- (Try different ML model)

