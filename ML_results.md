
# Results ML

The rationale behind using machine learning methods is to build a model that predicts PSD or not based on the gene expression profile of a patient. Training the model with an L1-penalty (penalize the use of more genes) induces that only "relevant" genes be included in the final model.

In all following models, 70% of the samples were used for training and the remaining 30% were used to validate the model.

The second study contains two datasets: one was used to train the model, and one to validate the model. Unfortunately, different platforms were used to acquire the gene expression profiles, leading to a different nomenclature of the gene names. As such, the second part of the second dataset could not be used in any analysis.

## Logistic Regression (scikit-learn LogisticRegression)

### Results
Training on only the second dataset led to a model averaging 89% classification accuracy using only 40 genes. This model however performed poorly on the first dataset, averaging only 56% classification accuracy, and thus averages 66% classification accuracy on the combined datasets. This model clearly overfits the second dataset and is of little value for further investigations.

Training on the combined datasets led to a model averaging 84% classification accuracy using only 50 genes. Training and validating on the first dataset only leads to a 73% classification accuracy, and 82% when only using the second dataset.

These 50 genes are:

| Name | 
| ------------- |
|1553243_at|
|1554675_a_at|
|1555097_a_at|
|1555591_at|
|1556209_at|
|1556617_a_at|
|1556706_at|
|1556786_at|
|1557161_at|
|1557890_at|
|1558857_at|
|1561003_at|
|1561365_at|
|1561578_s_at|
|1564010_at|
|1569975_at|
|202198_s_at|
|202933_s_at|
|205741_s_at|
|206835_at|
|207681_at|
|208526_at|
|209710_at|
|209997_x_at|
|212836_at|
|215636_at|
|219196_at|
|219758_at|
|220183_s_at|
|220783_at|
|221827_at|
|222324_at|
|223187_s_at|
|226137_at|
|228559_at|
|230514_s_at|
|231446_at|
|231487_at|
|233379_at|
|234110_at|
|235104_at|
|237597_at|
|238314_x_at|
|239951_at|
|240125_at|
|240750_at|
|244654_at|
|34408_at|
|56919_at|
|65635_at|

Further analysis: comparison with statistical results

### Implementation details
- 4-fold cross-validation using bootstrapping led to an optimal regularization constant C of 0.1 for L1-regularization.
- Maximum 100 iterations to train the model


## Stochastic Gradient Descent (scikit-learn SGDClassifier)

### Results
Training on only one dataset proved difficult, as the model selected huge numbers of relevant genes (around 30'000) to compensate for the limited amount of samples.

Training on the combined datasets however led to a model averaging 75% classification accuracy using only 39 genes. Training and validating on the first dataset only led to a 75% classification accuracy, and 77% when only using the second dataset.

Training a logistic regression model using only these 39 selected genes improved the average classification accuracy to 78%.

These 39 genes are:

| Name | 
| ------------- |
|1555707_at|
|1555827_at|
|1558920_at|
|1559470_at|
|1561778_at|
|1562933_at|
|1564052_at|
|1567702_at|
|1567860_at|
|201056_at|
|206702_at|
|206835_at|
|207681_at|
|209120_at|
|209710_at|
|210387_at|
|211298_s_at|
|213850_s_at|
|214882_s_at|
|216566_at|
|217621_at|
|219077_s_at|
|220783_at|
|223580_at|
|224396_s_at|
|231842_at|
|232752_at|
|238314_x_at|
|239186_at|
|239871_at|
|239975_at|
|240125_at|
|240488_at|
|240946_at|
|241001_at|
|242662_at|
|243283_at|
|243720_at|
|244505_at|

### Implementation details
- 4-fold cross-validation using bootstrapping led to an optimal regularization constant C of 0.1 for L1-regularization
- Hinge loss and elastic net penalty, with a L1-to-L2-ratio of 0.9
- Maximum 100 iterations of stochastic gradient descent

## Random Forest (scikit-learn RandomForestClassifier)

### Results
Random forests do not automatically select features in the same sense as the previous models. If the trees used are deep enough, a random forest will eventually use all features. 

Restricting the depth of the trees, an average of 85% classification accuracy was achieved on the combined datasets using 3949 genes. The number of used genes is too large and thus further analyses using random forests were abandoned.

### Implementation details
- Maximal depth of 20 nodes for each tree
- 100 estimators (trees) in each forest

## Technicalities
- For any approach, 100 models of the same type were trained, and the best one was saved if it used less than 100 genes. This process was repeated multiple times.
- On single dataset models, only genes with expression values of more than 150 in all samples were used to build the model
- On combined dataset models, all genes were used to build the model

# Analysis

## Logistic Regression vs Stochastic Gradient Descent

6 genes were found to be relevant in both models. These genes are

| Matching genes| 
| ------------- |
|206835_at|
|207681_at|
|209710_at <--|
|220783_at|
|238314_x_at|
|240125_at|   

## TODO
- Transfer to regular gene name
- Compare with statistical analysis
- Compare with ongoing research
- Statement about improving the paper's model
- No age or batch number as feature

