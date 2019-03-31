Geneset Enrichment Analysis
================
Virginia Pichler
March 30, 2019

#### Input metadata and expression data cleaned in the [pre\_processing\_data.Rmd](https://github.com/STAT540-UBC/Repo_team_Y0ung-parents_W2019/blob/master/code/Pre_processing.rmd)

#### Load and format data

``` r
#load in feature data
geo_GSE25507 <- getGEO("GSE25507", GSEMatrix = TRUE)
```

    ## Found 1 file(s)

    ## GSE25507_series_matrix.txt.gz

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   ID_REF = col_character()
    ## )

    ## See spec(...) for full column specifications.

    ## File stored at:

    ## C:\Users\Ginny\AppData\Local\Temp\RtmpYXSEW8/GPL570.soft

    ## Warning: 62 parsing failures.
    ##   row     col           expected    actual         file
    ## 54614 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## 54615 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## 54616 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## 54617 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## 54618 SPOT_ID 1/0/T/F/TRUE/FALSE --Control literal data
    ## ..... ....... .................. ......... ............
    ## See problems(...) for more details.

``` r
new_geo_GSE25507<- geo_GSE25507[[1]]

# pull gene symbol from feature data
symbol<-fData(new_geo_GSE25507)
symbol<-symbol[, c("ID", "Gene Symbol")]
symbol<-symbol%>%arrange(ID)


# load in expression data and metadata

Meta_data<-readRDS('Meta_data.rds')
combine_norm<-readRDS('combine_norm.rds')

# add column for sample id in metadata
Meta_data$sample_id <- rownames(Meta_data)

# move sample id column to first position in metadata
Meta_data <- Meta_data %>% select(sample_id, everything())
rownames(Meta_data) <- c()
```

#### Create linear model

``` r
# check equivalence of samples in metadata and expression data
Meta_data$sample_id==names(combine_norm)
```

    ## logical(0)

``` r
#design matrix
Meta_data_design <-Meta_data[,c("sample_id", "diagnosis" )]
designMatrix <- model.matrix(~diagnosis, Meta_data_design)


# lmfit 
lmfit <- lmFit(combine_norm, designMatrix)

# run ebayes to calculate moderated t-statistics
lmfit_ebayes <- eBayes(lmfit)
```

#### Select, format and annotate genes for enrichment analysis

``` r
# select all genes and reformat and arrange data frame
allGenes <- topTable(lmfit_ebayes, number = Inf)
```

    ## Removing intercept from test coefficients

``` r
allGenes$ID <- rownames(allGenes)
allGenes<-allGenes%>%arrange(ID)

# annotate with gene symbols
allGenes<-full_join(allGenes, symbol, by = "ID")
colnames(allGenes)[8] <- "gene"
allGenes<-unique(allGenes)
allGenes<-na.omit(allGenes)

# subset data frame to gene and logFC values.
subsettedallGenes<-allGenes[, c("gene", "logFC")]
subsettedallGenes<-unique(subsettedallGenes)
```

#### Geneset Enrichment Analysis

``` r
# retrieve gene multifunctionality scores
urlfile<-'https://raw.githubusercontent.com/STAT540-UBC/STAT540-UBC.github.io/master/seminars/seminars_winter_2019/seminar10/data/gene_multifunctionality_scores.csv'

gene_multifunctionality_scores<-read.csv(urlfile)


# check if gene list has multifunctional bias
mergedData <- subsettedallGenes %>% inner_join(gene_multifunctionality_scores, by = "gene")

rankMfCor <- cor(abs(mergedData$logFC), mergedData$MF.score, method = "spearman")


# plot Spearman's correlation 
mergedData %>%
  ggplot(aes(x = abs(logFC), y = log10(MF.score))) + 
  geom_point(alpha = 0.1) +
  geom_smooth(method = "lm", se = FALSE) +
  ggtitle(paste0("r = ", rankMfCor))
```

![](Geneset_Enrichment_Analysis_files/figure-markdown_github/unnamed-chunk-4-1.png)

``` r
### geneset enrichment analysis

#retrieve GO.xml file from ermineR
if (!file.exists("GO.xml")) { goToday("GO.xml") }


# convert scores to absolute values
ermineInputGeneScores <- subsettedallGenes %>% 
  mutate(absolute_logFC = abs(logFC)) %>% 
  select(gene, absolute_logFC) %>% 
  na.omit() %>% 
  as.data.frame() %>% 
  arrange(desc(absolute_logFC))

rownames(ermineInputGeneScores) <- make.names(subsettedallGenes[,1], unique = TRUE)

# enrichment with Precision-Recall method
enrichmentResult <- precRecall(scores = ermineInputGeneScores, 
                               scoreColumn = 2, # scores 
                               bigIsBetter = TRUE, # rank large logFC higher
                               annotation = "Generic_human", # ermineR Generic_human annotation file
                               aspects = "B", # biological processes 
                               iterations = 10000, # 10K sampling iterations so results are stable
                               geneSetDescription = "GO.xml") # GO XML file


# arrange results by MFP value
enrichmentResult$results %>% arrange(MFPvalue)
```

    ## # A tibble: 3,494 x 12
    ##    Name  ID    NumProbes NumGenes RawScore    Pval CorrectedPvalue MFPvalue
    ##    <chr> <chr>     <dbl>    <dbl>    <dbl>   <dbl>           <dbl>    <dbl>
    ##  1 posi~ GO:2~        62       62   0.0225 0.0001            0.335   0.0001
    ##  2 nucl~ GO:0~        36       36   0.0308 0.0004            0.670   0.0002
    ##  3 erro~ GO:0~        20       20   0.0517 0.0008            0.893   0.0008
    ##  4 tran~ GO:0~        72       72   0.0197 0.0008            0.670   0.0008
    ##  5 DNA ~ GO:0~       191      191   0.0194 0.00120           0.804   0.001 
    ##  6 regu~ GO:2~       102      102   0.0178 0.0013            0.726   0.0011
    ##  7 thyr~ GO:0~        21       21   0.0458 0.0014            0.670   0.0014
    ##  8 telo~ GO:0~        22       22   0.0473 0.0014            0.586   0.0014
    ##  9 nucl~ GO:0~        36       36   0.0302 0.0017            0.633   0.0017
    ## 10 telo~ GO:0~        93       93   0.0170 0.0027            0.476   0.0017
    ## # ... with 3,484 more rows, and 4 more variables: CorrectedMFPvalue <dbl>,
    ## #   Multifunctionality <dbl>, `Same as` <chr>, GeneMembers <chr>

``` r
# scatterplot of multifunctionality adjustment to GO terms
enrichmentResult$results %>% 
  ggplot(aes(x = -log10(Pval), y = -log10(MFPvalue))) +
  geom_point(alpha = 0.2)
```

![](Geneset_Enrichment_Analysis_files/figure-markdown_github/unnamed-chunk-4-2.png)

``` r
# top ten GO terms with largest adjustments 
enrichmentResult$results %>% 
  select(Name, ID, Pval, MFPvalue) %>% 
  mutate(neg_log_pvalue = -log10(Pval),
         neg_log_mfpvalue = -log10(MFPvalue)) %>% 
  mutate(log_pvalue_change = neg_log_mfpvalue - neg_log_pvalue) %>% 
  arrange(desc(abs(log_pvalue_change))) %>% 
  head(10) %>% 
  kable()
```

| Name                                                         | ID           |    Pval|  MFPvalue|  neg\_log\_pvalue|  neg\_log\_mfpvalue|  log\_pvalue\_change|
|:-------------------------------------------------------------|:-------------|-------:|---------:|-----------------:|-------------------:|--------------------:|
| positive regulation of immune effector process               | <GO:0002699> |  0.2156|    0.0790|         0.6663512|           1.1023729|            0.4360217|
| response to corticosteroid                                   | <GO:0031960> |  0.0903|    0.0344|         1.0443122|           1.4634416|            0.4191293|
| regulation of blood vessel size                              | <GO:0050880> |  0.1961|    0.0774|         0.7075224|           1.1112590|            0.4037366|
| regulation of muscle tissue development                      | <GO:1901861> |  0.2470|    0.0982|         0.6073030|           1.0078885|            0.4005855|
| regulation of carbohydrate metabolic process                 | <GO:0006109> |  0.1272|    0.0509|         0.8955129|           1.2932822|            0.3977693|
| regulation of muscle contraction                             | <GO:0006937> |  0.1272|    0.0509|         0.8955129|           1.2932822|            0.3977693|
| vascular process in circulatory system                       | <GO:0003018> |  0.1821|    0.0734|         0.7396901|           1.1343039|            0.3946139|
| response to ketone                                           | <GO:1901654> |  0.0745|    0.0301|         1.1278437|           1.5214335|            0.3935898|
| regulation of cytokine secretion                             | <GO:0050707> |  0.2912|    0.1177|         0.5358086|           0.9292235|            0.3934149|
| negative regulation of establishment of protein localization | <GO:1904950> |  0.2339|    0.0948|         0.6309698|           1.0231917|            0.3922219|
