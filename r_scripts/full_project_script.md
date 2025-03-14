# 1. Goal of the predictor

Accurate estimation of prognosis is essential for clinical decision
making and care for breast cancer patients. Patients with same clinical
characteristics (e.g., tumor stage) can respond to treatments
differently and experience different survival outcomes, which might be
attributable to complex variation in individual gene expression. Hormone
and radiation therapy are common treatments for breast cancer and there
are many patients who receive both. Recent studies raised the concern
that hormonal therapy may reduce the efficacy of radiation [(Cecchini et
al.,
2015)](https://www.ncbi.nlm.nih.govpmc/articles/PMC4659580/#:~:text=However%2C%20when%20comparing%20radiation%20alone,outcomes%2C%20including%20survival%20and%20recurrence.),
and it is not clear whether accurate estimation of survival can be
achieved depending on patient treatment assignments. The **goal** of the
current analysis is to create predictors for survival among patients who
received only hormone, only radiation, or both hormone and radiation
therapy.

# 2. Study design

The data is from the Molecular Taxonomy of Breast Cancer International
Consortium (METABRIC) database and accessed from [Kaggle (clink for
link)](https://www.kaggle.com/datasets/raghadalharbi/breast-cancer-gene-expression-profiles-metabric).
The prospective data consists of 31 clinical features and 506 genetic
attributes collected from almost 2000 breast cancer patients. These
patients were followed up and the time (months) to either censoring or
death were collected.

A part of the genetic attributes are numeric z-scores representing the
number of standard deviations in mRNA expression from the mean
expression in a reference population, which measures whether a gene is
increased or decreased relative to tumor-free people or patients of
other tumors. The rest of the genetic features are binary indicators of
mutation, which is very sparse. The clinical attributes include 31
features such as cancer stage, treatment, tumor size, etc.

The structure of the data lead me to conduct a survival analysis. I also
decided to drop the mutation features in the data for they are too
sparse. One of the clinical variable summarized the number of mutations
of each patient, so I felt this is a good enough representation of the
mutation data.

# 3. Algorithms

Because I have time-to-event data, the algorithms that I tried are the
**Cox proportional hazard regression**, **Cox proportional hazard
regression with regularization (lasso, ridge, alpha=0.5)**, and **random
survival forest**. I choose the C-index as my loss function for its good
interpretability like AUROC.

For the regularized cox models, I plan to use cross-validation to tune
the lambda value and use the lambda value that gives minimum mean
cross-validated error. For the random survival forest model, I use the
parameters that were suggested by [Pickett and colleagues
(2021)](https://bmcmedresmethodol.biomedcentral.com/articles/10.1186/s12874-021-01375-x#Tab1).

``` r
#original dataset dimension
df<- read.csv("raw_data/METABRIC_RNA_Mutation 2.csv") 
dim(df)
```

    ## [1] 1904  693

# 4.Feature engenieering and exploratory data analysis

## 4.1 Data dimension & key variable distributions

The orginal dataset contains 1904 observations and 693 features. Among
the features, there are 31 clinical variables. There are 331 mRNA level
z-scores (numeric) for 331 genes and mutation (binary) for 175 genes.
Data on mutation of genes are very sparse, so that I decided to just
focus on clinical features and mRNA data.

My exploration indicates that 801 deaths occurred during the follow up
period. The average age of the cohort if 61. Patients in this cohort
receive three types of therapies and combinations of them. The data also
shows that very few patients have mutation on more than 15 genes. I also
found that **missing values** were only partially coded as NA. Some of
them are empty cells.

``` r
#binary overall survival status
table(df$overall_survival)
```

    ## 
    ##    0    1 
    ## 1103  801

``` r
#age 
summary(df$age_at_diagnosis)
```

    ##    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
    ##   21.93   51.38   61.77   61.09   70.59   96.29

``` r
#cancer type 
table(df$cancer_type_detailed)
```

    ## 
    ##                                           
    ##                                        15 
    ##                                    Breast 
    ##                                        17 
    ##          Breast Invasive Ductal Carcinoma 
    ##                                      1500 
    ##         Breast Invasive Lobular Carcinoma 
    ##                                       142 
    ##  Breast Invasive Mixed Mucinous Carcinoma 
    ##                                        22 
    ## Breast Mixed Ductal and Lobular Carcinoma 
    ##                                       207 
    ##                 Metaplastic Breast Cancer 
    ##                                         1

``` r
#treatment 
table(df$chemotherapy)
```

    ## 
    ##    0    1 
    ## 1508  396

``` r
table(df$hormone_therapy)
```

    ## 
    ##    0    1 
    ##  730 1174

``` r
table(df$radio_therapy)
```

    ## 
    ##    0    1 
    ##  767 1137

``` r
table(df$chemotherapy, df$hormone_therapy, df$radio_therapy)
```

    ## , ,  = 0
    ## 
    ##    
    ##       0   1
    ##   0 289 405
    ##   1  45  28
    ## 
    ## , ,  = 1
    ## 
    ##    
    ##       0   1
    ##   0 228 586
    ##   1 168 155

``` r
#number of mutation 
table(df$mutation_count)
```

    ## 
    ##   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20 
    ## 107 193 229 248 268 232 168 121  90  61  38  25  17  16  11   8   5   2   2   4 
    ##  21  22  23  24  26  28  30  40  46  80 
    ##   1   4   2   1   1   1   1   1   1   1

``` r
#missing values ranked; not all are correctly represented
missing<- data.frame(sapply(df, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df..function.x..sum.is.na.x...))  
head(missing, 10)
```

    ##                           sapply.df..function.x..sum.is.na.x...
    ## tumor_stage                                                 501
    ## neoplasm_histologic_grade                                    72
    ## mutation_count                                               45
    ## tumor_size                                                   20
    ## patient_id                                                    0
    ## age_at_diagnosis                                              0
    ## type_of_breast_surgery                                        0
    ## cancer_type                                                   0
    ## cancer_type_detailed                                          0
    ## cellularity                                                   0

# 5. Data pre-processing

In my data pre-processing, I follow the steps below:

-   Explore missing values and correctly code missing as NAs
-   Manually drop features that will not be useful for this task
-   Multiple imputation by chained equation
-   Create categorical variables to stratify patients based on treatment

## 5.1 Explore missing values and replace missing with NAs

``` r
#there isn't missing values in the outcome 
table(is.na(df$overall_survival))
```

    ## 
    ## FALSE 
    ##  1904

``` r
table(is.na(df$overall_survival_months))
```

    ## 
    ## FALSE 
    ##  1904

``` r
#replace empty cells with NA 
df<- df %>% 
  replace_with_na_all(condition= ~.x=="")
#a subset of clinical variables for viewing purpose
clinical<- clinical %>% 
  replace_with_na_all(condition= ~.x=="")
```

## 5.2 Drop clinical features that will not be useful

``` r
#clinical feature list
str(clinical)
```

    ## tibble [1,904 × 31] (S3: tbl_df/tbl/data.frame)
    ##  $ patient_id                    : int [1:1904] 0 2 5 6 8 10 14 22 28 35 ...
    ##  $ age_at_diagnosis              : num [1:1904] 75.7 43.2 48.9 47.7 77 ...
    ##  $ type_of_breast_surgery        : chr [1:1904] "MASTECTOMY" "BREAST CONSERVING" "MASTECTOMY" "MASTECTOMY" ...
    ##  $ cancer_type                   : chr [1:1904] "Breast Cancer" "Breast Cancer" "Breast Cancer" "Breast Cancer" ...
    ##  $ cancer_type_detailed          : chr [1:1904] "Breast Invasive Ductal Carcinoma" "Breast Invasive Ductal Carcinoma" "Breast Invasive Ductal Carcinoma" "Breast Mixed Ductal and Lobular Carcinoma" ...
    ##  $ cellularity                   : chr [1:1904] NA "High" "High" "Moderate" ...
    ##  $ chemotherapy                  : int [1:1904] 0 0 1 1 1 0 1 0 0 0 ...
    ##  $ pam50_._claudin.low_subtype   : chr [1:1904] "claudin-low" "LumA" "LumB" "LumB" ...
    ##  $ cohort                        : num [1:1904] 1 1 1 1 1 1 1 1 1 1 ...
    ##  $ er_status_measured_by_ihc     : chr [1:1904] "Positve" "Positve" "Positve" "Positve" ...
    ##  $ er_status                     : chr [1:1904] "Positive" "Positive" "Positive" "Positive" ...
    ##  $ neoplasm_histologic_grade     : num [1:1904] 3 3 2 2 3 3 2 2 3 2 ...
    ##  $ her2_status_measured_by_snp6  : chr [1:1904] "NEUTRAL" "NEUTRAL" "NEUTRAL" "NEUTRAL" ...
    ##  $ her2_status                   : chr [1:1904] "Negative" "Negative" "Negative" "Negative" ...
    ##  $ tumor_other_histologic_subtype: chr [1:1904] "Ductal/NST" "Ductal/NST" "Ductal/NST" "Mixed" ...
    ##  $ hormone_therapy               : int [1:1904] 1 1 1 1 1 1 1 1 1 0 ...
    ##  $ inferred_menopausal_state     : chr [1:1904] "Post" "Pre" "Pre" "Pre" ...
    ##  $ integrative_cluster           : chr [1:1904] "4ER+" "4ER+" "3" "9" ...
    ##  $ primary_tumor_laterality      : chr [1:1904] "Right" "Right" "Right" "Right" ...
    ##  $ lymph_nodes_examined_positive : num [1:1904] 10 0 1 3 8 0 1 1 1 0 ...
    ##  $ mutation_count                : num [1:1904] NA 2 2 1 2 4 4 1 4 5 ...
    ##  $ nottingham_prognostic_index   : num [1:1904] 6.04 4.02 4.03 4.05 6.08 ...
    ##  $ oncotree_code                 : chr [1:1904] "IDC" "IDC" "IDC" "MDLC" ...
    ##  $ overall_survival_months       : num [1:1904] 140.5 84.6 163.7 164.9 41.4 ...
    ##  $ overall_survival              : int [1:1904] 1 1 0 1 0 0 1 0 0 0 ...
    ##  $ pr_status                     : chr [1:1904] "Negative" "Positive" "Positive" "Positive" ...
    ##  $ radio_therapy                 : int [1:1904] 1 1 0 1 1 1 1 1 1 0 ...
    ##  $ X3.gene_classifier_subtype    : chr [1:1904] "ER-/HER2-" "ER+/HER2- High Prolif" NA NA ...
    ##  $ tumor_size                    : num [1:1904] 22 10 15 25 40 31 10 29 16 28 ...
    ##  $ tumor_stage                   : num [1:1904] 2 1 2 2 2 4 2 2 2 2 ...
    ##  $ death_from_cancer             : chr [1:1904] "Living" "Living" "Died of Disease" "Living" ...

``` r
#drop redundant, not useful, outcome features.. 
drop<-c("cancer_type", "cohort","patient_id", "er_status_measured_by_ihc", "her2_status_measured_by_snp6","nottingham_prognostic_index","death_from_cancer", "overall_survival", "overall_survival_months")
outcome<- df[, colnames(df) %in% c("overall_survival", "overall_survival_months")]
df_select<-df[, !colnames(df) %in% drop]

#drop mutations in genes because too sparse 
df_select<- df_select %>%
  select(-ends_with("_mut"))

#dims after dropping 
dim(df_select)
```

    ## [1] 1904  511

``` r
dim(outcome)
```

    ## [1] 1904    2

## 5.3 Multiple imupatation

Now that all missing are correctly coded, I first rank the features by
numbers of missingness. I can see that the missing values are all on the
clinical features. I use multiple imputation with chained equations
depending on the type of the features. The data table below looks like
there is no missingness.

The imputation process is muted here because I saved the imputed dataset
to save time. But the code is shown here.

``` r
#rank by missing #
missing<- data.frame(sapply(df_select, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df_select..function.x..sum.is.na.x...)) 
head(missing, 15)
```

    ##                                sapply.df_select..function.x..sum.is.na.x...
    ## tumor_stage                                                             501
    ## X3.gene_classifier_subtype                                              204
    ## primary_tumor_laterality                                                106
    ## neoplasm_histologic_grade                                                72
    ## cellularity                                                              54
    ## mutation_count                                                           45
    ## type_of_breast_surgery                                                   22
    ## tumor_size                                                               20
    ## cancer_type_detailed                                                     15
    ## tumor_other_histologic_subtype                                           15
    ## oncotree_code                                                            15
    ## age_at_diagnosis                                                          0
    ## chemotherapy                                                              0
    ## pam50_._claudin.low_subtype                                               0
    ## er_status                                                                 0

``` r
#convert all character features to factor for imputation
df_select[sapply(df_select, is.character)] <- lapply(df_select[sapply(df_select, is.character)], 
                                                           as.factor)
#imputation (excluding the outcome features)
set.seed(777)
imp<- mice(df_select[, !colnames(df_select) %in% c("overall_survival_months","overall_survival")], m=3, method = "pmm")
df_select_imp<- complete(imp)

#add outcome back
df_select_imp<- data.frame(df_select_imp, outcome)
```

``` r
#save this imputed data so that I don't have to run again 
save(df_select_imp, file="/processed_data/df_select_imp.rda")
```

``` r
load("./processed_data/df_select_imp.rda")
#check missing after imputation (ranked)
missing<- data.frame(sapply(df_select_imp, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df_select_imp..function.x..sum.is.na.x...))
head(missing, 10)
```

    ##                                sapply.df_select_imp..function.x..sum.is.na.x...
    ## age_at_diagnosis                                                              0
    ## type_of_breast_surgery                                                        0
    ## cancer_type_detailed                                                          0
    ## cellularity                                                                   0
    ## chemotherapy                                                                  0
    ## pam50_._claudin.low_subtype                                                   0
    ## er_status                                                                     0
    ## neoplasm_histologic_grade                                                     0
    ## her2_status                                                                   0
    ## tumor_other_histologic_subtype                                                0

## 5.4 Feature selection based on correlation

I dropped one of a pair of features if the pair correlation is over 0.7.
I ended up with 497 features now.

``` r
#correlation matrix of numeric features 
num_cols<-select_if(df_select_imp, is.numeric)
cor<- melt(round(x= cor(num_cols), digits = 2)) %>% 
  arrange(desc(abs(value))) %>% 
  filter(Var1!=Var2)
 
#remove |r|>0.7 
big_cor<- cor %>% 
  filter(abs(value)>0.7)
drop_r<- as.vector(big_cor[,1]) %>% 
  unique()
#drop these features 
df_select_imp<- df_select_imp[, !colnames(df_select_imp)%in% drop_r]
dim(df_select_imp)
```

    ## [1] 1904  497

## 5.5 Stratify based on treatment

There three main treatment types in this data, which are chemo,
radiation, and hormone. Many patients have more than 1 treatment.
Corresponding to my stated purpose, I decide to limit my data to
patients who receive hormone and/or radiation therapy. I ended up with
**1219 patients and 495** features. These patients are further
stratified into hormone only, radiation only, and both for the purpose
of creating different predictors.

``` r
# chemotherapy, hormone_therapy, radio_therapy
df_select_imp<- df_select_imp %>% 
  mutate(treatment= case_when(chemotherapy==1&hormone_therapy==1&radio_therapy==1 ~"all",
                                 chemotherapy==1&hormone_therapy==1&radio_therapy==0 ~"che_hor",
                                 chemotherapy==1&hormone_therapy==0&radio_therapy==1 ~"che_rad", 
                                 chemotherapy==1&hormone_therapy==0&radio_therapy==0 ~"che", 
                                 chemotherapy==0&hormone_therapy==1&radio_therapy==1 ~"hor_rad", 
                                 chemotherapy==0&hormone_therapy==0&radio_therapy==1 ~"rad", 
                                 chemotherapy==0&hormone_therapy==1&radio_therapy==0 ~ "hor", 
                                 chemotherapy==0&hormone_therapy==0&radio_therapy==0 ~"none"
                                 )) %>% 
  filter(treatment %in% c("hor", "hor_rad", "rad")) %>% 
  mutate(treatment= factor(treatment, levels= c("hor_rad", "hor", "rad"))) %>% 
  select(-chemotherapy, -hormone_therapy, -radio_therapy) %>% 
  rename(time= overall_survival_months, status= overall_survival)

dim(df_select_imp)
```

    ## [1] 1219  495

``` r
table(df_select_imp$treatment)
```

    ## 
    ## hor_rad     hor     rad 
    ##     586     405     228

``` r
#drop "chemotherapy", "hormone_therapy", "radio_therapy")) 
```

# 6. Model training and evaluation

I generally followed the following steps:

-   Explore survival curves of patients who had different treatment
-   Fit different models using the afromentioned algorithms and
    parameters tuning for patients of each of the 3 types of treatment
-   Evaluate the 10-fold cross-validated C-index
-   Repeat the last two steps using the full data set of all patients
    who have had hormone and/or radiation therapy
-   Compare performance

## 6.1 Survival curves of patients who received different treatment

We can see that the unadjusted survival curves are significantly
different among patients who receive different treatments. Echoing the
concern about reduced efficacy of using both hormone and radiation
therapy, we do see the worst survival for patients who receive both
treatment. But this is likely confounded by other factors such as tumor
stage and type.

``` r
cox1<- survfit(Surv(time, status)~ treatment, data=df_select_imp)
ggsurvplot(cox1, data= df_select_imp,
           title= "Survival curves of breast cancer patients by treatment group", 
           legend.title="",
           xlab= "Months since study entry",
           legend.labs=c("Hormone + Radiation Therapy", "Only Hormone Therapy", "Only Radation Therapy"),
           pval=T, conf.int = T, risk.table = T, tables.height = 0.2, tables.theme = theme_cleantable(), ggtheme=theme_bw(),
          )
```

![](analysis0601_files/figure-markdown_github/unnamed-chunk-14-1.png)

``` r
#subset data 
hor_rad<- df_select_imp %>% 
  filter(treatment=="hor_rad") %>% 
  select(-treatment)
hor<- df_select_imp %>% 
  filter(treatment=="hor") %>% 
  select(-treatment)
rad<- df_select_imp %>% 
  filter(treatment=="rad") %>% 
  select(-treatment)
```

## 6.2 Predictors for patients who received hormone and radiation therapy

``` r
set.seed(777)
N<- nrow(hor_rad)
V<- 10
folds<- split(1:N, rep(1:10, length=N))
eval1<- matrix(NA, nrow=V, ncol = 5)
colnames(eval1)<-c("coxph", "lasso", "ridge", "el", "rf")

#coxph
for (v in 1:V){
  train<- hor_rad[-folds[[v]], ]
  test<- hor_rad[folds[[v]], ]
  fit_cv_cph<- coxph(Surv(time, status)~., data= train)
  pred_cph<- predict(fit_cv_cph, 
                     type= "lp", 
                     newdata= test)
  eval1[v,1]<- concordance(Surv(time, status)~ pred_cph, data= test, reverse = TRUE)$concordance
}
```

``` r
#lasso
Y<- data.matrix(hor_rad[c("time", "status")])
X<- data.matrix(hor_rad[, !colnames(hor_rad)%in% c("time", "status")])

for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=1)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval1[v,2]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#ridge
for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=0)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval1[v,3]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#el
for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=0.5)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval1[v,4]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}
```

``` r
#rf
for (v in 1:V){
  train<- hor_rad[-folds[[v]], ]
  test<- hor_rad[folds[[v]], ]
  fit_cv_rf<- rfsrc(Surv(time, status)~., data= train, 
                    ntree=500, 
                    samptype = "swr")#default for other pars, see reference 
  eval1[v,5]<- 1- tail(predict.rfsrc(fit_cv_rf, newdata = test)$err.rate,1)
}
```

## 6.3 Predictors for patients who received only hormone therapy

``` r
set.seed(555)

N<- nrow(hor)
V<- 10
folds<- split(1:N, rep(1:V, length=N))
eval2<- matrix(NA, nrow=V, ncol = 5)
colnames(eval2)<-c("coxph", "lasso", "ridge", "el", "rf")

#coxph
for (v in 1:V){
  train<- hor[-folds[[v]], ]
  test<- hor[folds[[v]], ]
  fit_cv_cph<- coxph(Surv(time, status)~., data= train)
  pred_cph<- predict(fit_cv_cph, 
                     type= "lp", 
                     newdata= test)
  eval2[v,1]<- concordance(Surv(time, status)~ pred_cph, data= test, reverse = TRUE)$concordance
}
```

``` r
#lasso
Y<- data.matrix(hor[c("time", "status")])
X<- data.matrix(hor[, !colnames(hor)%in% c("time", "status")])

for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=1)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval2[v,2]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#ridge
for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=0)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval2[v,3]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#el
for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=0.5)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval2[v,4]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}
```

``` r
#rf
for (v in 1:V){
  train<- hor[-folds[[v]], ]
  test<- hor[folds[[v]], ]
  fit_cv_rf<- rfsrc(Surv(time, status)~., data= train, 
                    ntree=500, 
                    samptype = "swr")#default for other pars, see reference 
  eval2[v,5]<- 1- tail(predict.rfsrc(fit_cv_rf, newdata = test)$err.rate,1)
}
```

## 6.4 Predictors for patients who received only radiation therapy

``` r
set.seed(333)
N<- nrow(rad)
V<- 10
folds<- split(1:N, rep(1:V, length=N))
eval3<- matrix(NA, nrow=V, ncol = 5)
colnames(eval3)<-c("coxph", "lasso", "ridge", "el", "rf")

#coxph
for (v in 1:V){
  train<- rad[-folds[[v]], ]
  test<- rad[folds[[v]], ]
  fit_cv_cph<- coxph(Surv(time, status)~., data= train)
  pred_cph<- predict(fit_cv_cph, 
                     type= "lp", 
                     newdata= test)
  eval3[v,1]<- concordance(Surv(time, status)~ pred_cph, data= test, reverse = TRUE)$concordance
}
```

``` r
#lasso
Y<- data.matrix(rad[c("time", "status")])
X<- data.matrix(rad[, !colnames(rad)%in% c("time", "status")])

for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=1)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval3[v,2]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#ridge

for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=0)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval3[v,3]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#el

for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=0.5)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval3[v,4]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#rf
for (v in 1:V){
  train<- rad[-folds[[v]], ]
  test<- rad[folds[[v]], ]
  fit_cv_rf<- rfsrc(Surv(time, status)~., data= train, 
                    ntree=500, 
                    samptype = "swr")#default for other pars, see reference 
  eval3[v,5]<- 1- tail(predict.rfsrc(fit_cv_rf, newdata = test)$err.rate,1)
}
```

## 6.5 Predictors for all patients

``` r
set.seed(111)
N<- nrow(df_select_imp)
V<- 10
folds<- split(1:N, rep(1:V, length=N))
eval_all<- matrix(NA, nrow=V, ncol = 5)
colnames(eval_all)<-c("coxph", "lasso", "ridge", "el", "rf")

#coxph
for (v in 1:V){
  train<- df_select_imp[-folds[[v]], ]
  test<- df_select_imp[folds[[v]], ]
  fit_cv_cph<- coxph(Surv(time, status)~., data= train)
  pred_cph<- predict(fit_cv_cph, 
                     type= "lp", 
                     newdata= test)
  eval_all[v,1]<- concordance(Surv(time, status)~ pred_cph, data= test, reverse = TRUE)$concordance
}
```

``` r
Y<- data.matrix(df_select_imp[c("time", "status")])
X<- data.matrix(df_select_imp[, !colnames(df_select_imp)%in% c("time", "status")])
#lasso
for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=1)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval_all[v,2]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#ridge
for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=0)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval_all[v,3]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}

#el
for (v in 1:V){
  train_x<- X[-folds[[v]], ]
  test_x<- X[folds[[v]], ]
  train_y<- Y[-folds[[v]], ]
  test_y<- Y[folds[[v]], ]
  fit_cv_lasso<- cv.glmnet(y= train_y,
                           x= train_x, 
                           family= "cox", 
                           alpha=0.5)
  pred_lasso<- predict(fit_cv_lasso$glmnet.fit, 
                       newx= test_x, 
                       s= fit_cv_lasso$lambda.min, 
                       type="link")
  df_test<- data.frame(test_y, test_x)
  eval_all[v,4]<- concordance(Surv(time, status)~ pred_lasso, 
                          data= df_test, reverse = TRUE)$concordance
}
```

``` r
#rf
for (v in 1:V){
  train<- df_select_imp[-folds[[v]], ]
  test<- df_select_imp[folds[[v]], ]
  fit_cv_rf<- rfsrc(Surv(time, status)~., data= train, 
                    ntree=500, 
                    samptype = "swr")#default for other pars, see reference 
  eval_all[v,5]<- 1- tail(predict.rfsrc(fit_cv_rf, newdata = test)$err.rate,1)
}
```

## 6.6 Summary of cross-validated performance (C-index)

Based on the results, ridge predictor had the best performance in the
predictor for patients who had hormone and radiation therapy. Lasso
predictor had the best performance in the predictors for patients who
only had hormone therapy. Random survival forest had the best
performance for the predictor for patients who only had radiation
therapy and overall patients.

The overall performance is not that great, which might related to the
fact that there are too many features vs. patients number. The
performance using regularized regression and random survival forest are
also close to each other. They all outperform the traditional cox
proportional regression.

``` r
#patients had both hormone and radiation therapy
summary(eval1)
```

    ##      coxph            lasso            ridge              el        
    ##  Min.   :0.3642   Min.   :0.5775   Min.   :0.6285   Min.   :0.5855  
    ##  1st Qu.:0.4544   1st Qu.:0.5962   1st Qu.:0.6502   1st Qu.:0.6059  
    ##  Median :0.4879   Median :0.6663   Median :0.6669   Median :0.6678  
    ##  Mean   :0.4904   Mean   :0.6576   Mean   :0.6888   Mean   :0.6629  
    ##  3rd Qu.:0.5349   3rd Qu.:0.7014   3rd Qu.:0.7314   3rd Qu.:0.6957  
    ##  Max.   :0.5785   Max.   :0.7692   Max.   :0.7713   Max.   :0.7817  
    ##        rf        
    ##  Min.   :0.5996  
    ##  1st Qu.:0.6486  
    ##  Median :0.6628  
    ##  Mean   :0.6661  
    ##  3rd Qu.:0.7014  
    ##  Max.   :0.7271

``` r
#patients only had hormone therapy
summary(eval2)
```

    ##      coxph            lasso            ridge              el        
    ##  Min.   :0.4246   Min.   :0.5474   Min.   :0.4772   Min.   :0.5404  
    ##  1st Qu.:0.4627   1st Qu.:0.6252   1st Qu.:0.5886   1st Qu.:0.6201  
    ##  Median :0.5571   Median :0.6493   Median :0.6187   Median :0.6402  
    ##  Mean   :0.5557   Mean   :0.6702   Mean   :0.6337   Mean   :0.6725  
    ##  3rd Qu.:0.6368   3rd Qu.:0.7218   3rd Qu.:0.6962   3rd Qu.:0.7246  
    ##  Max.   :0.7219   Max.   :0.8411   Max.   :0.8013   Max.   :0.8543  
    ##        rf        
    ##  Min.   :0.5000  
    ##  1st Qu.:0.5537  
    ##  Median :0.6542  
    ##  Mean   :0.6443  
    ##  3rd Qu.:0.7406  
    ##  Max.   :0.7748

``` r
#patients only had radiation therapy
summary(eval3)
```

    ##      coxph            lasso            ridge              el        
    ##  Min.   :0.3548   Min.   :0.4750   Min.   :0.4757   Min.   :0.4333  
    ##  1st Qu.:0.4111   1st Qu.:0.5479   1st Qu.:0.5371   1st Qu.:0.5352  
    ##  Median :0.4541   Median :0.5855   Median :0.5941   Median :0.5969  
    ##  Mean   :0.4767   Mean   :0.5766   Mean   :0.6023   Mean   :0.5778  
    ##  3rd Qu.:0.5152   3rd Qu.:0.6194   3rd Qu.:0.6473   3rd Qu.:0.6400  
    ##  Max.   :0.6600   Max.   :0.6723   Max.   :0.7692   Max.   :0.6723  
    ##        rf        
    ##  Min.   :0.4663  
    ##  1st Qu.:0.5191  
    ##  Median :0.6327  
    ##  Mean   :0.6166  
    ##  3rd Qu.:0.6618  
    ##  Max.   :0.7802

``` r
#among all patients
summary(eval_all)
```

    ##      coxph            lasso            ridge              el        
    ##  Min.   :0.5124   Min.   :0.5816   Min.   :0.5948   Min.   :0.5822  
    ##  1st Qu.:0.5334   1st Qu.:0.6420   1st Qu.:0.6559   1st Qu.:0.6444  
    ##  Median :0.5470   Median :0.6817   Median :0.6820   Median :0.6805  
    ##  Mean   :0.5648   Mean   :0.6775   Mean   :0.6862   Mean   :0.6776  
    ##  3rd Qu.:0.6016   3rd Qu.:0.7128   3rd Qu.:0.7115   3rd Qu.:0.7160  
    ##  Max.   :0.6543   Max.   :0.7483   Max.   :0.7615   Max.   :0.7476  
    ##        rf        
    ##  Min.   :0.5816  
    ##  1st Qu.:0.6658  
    ##  Median :0.6830  
    ##  Mean   :0.6867  
    ##  3rd Qu.:0.7262  
    ##  Max.   :0.7652

# 6.7 Feature importance plot of the random survival forest predictor among all patients

``` r
set.seed(500)
rf_1<- rfsrc(Surv(time, status)~., data= df_select_imp, ntree=500, 
                  samptype = "swr", importance=TRUE)
jk.obj<- subsample(rf_1)#default B=100 #of bootstrap
```

    ## [[32m[0m                                                  ]   1%[[32m█[0m                                                 ]   2%[[32m██[0m                                                ]   3%[[32m██[0m                                                ]   4%[[32m██[0m                                                ]   5%[[32m███[0m                                               ]   6%[[32m████[0m                                              ]   7%[[32m████[0m                                              ]   8%[[32m████[0m                                              ]   9%[[32m█████[0m                                             ]  10%[[32m██████[0m                                            ]  11%[[32m██████[0m                                            ]  12%[[32m██████[0m                                            ]  13%[[32m███████[0m                                           ]  14%[[32m████████[0m                                          ]  15%[[32m████████[0m                                          ]  16%[[32m████████[0m                                          ]  17%[[32m█████████[0m                                         ]  18%[[32m██████████[0m                                        ]  19%[[32m██████████[0m                                        ]  20%[[32m██████████[0m                                        ]  21%[[32m███████████[0m                                       ]  22%[[32m████████████[0m                                      ]  23%[[32m████████████[0m                                      ]  24%[[32m████████████[0m                                      ]  25%[[32m█████████████[0m                                     ]  26%[[32m██████████████[0m                                    ]  27%[[32m██████████████[0m                                    ]  28%[[32m██████████████[0m                                    ]  29%[[32m███████████████[0m                                   ]  30%[[32m████████████████[0m                                  ]  31%[[32m████████████████[0m                                  ]  32%[[32m████████████████[0m                                  ]  33%[[32m█████████████████[0m                                 ]  34%[[32m██████████████████[0m                                ]  35%[[32m██████████████████[0m                                ]  36%[[32m██████████████████[0m                                ]  37%[[32m███████████████████[0m                               ]  38%[[32m████████████████████[0m                              ]  39%[[32m████████████████████[0m                              ]  40%[[32m████████████████████[0m                              ]  41%[[32m█████████████████████[0m                             ]  42%[[32m██████████████████████[0m                            ]  43%[[32m██████████████████████[0m                            ]  44%[[32m██████████████████████[0m                            ]  45%[[32m███████████████████████[0m                           ]  46%[[32m████████████████████████[0m                          ]  47%[[32m████████████████████████[0m                          ]  48%[[32m████████████████████████[0m                          ]  49%[[32m█████████████████████████[0m                         ]  50%[[32m██████████████████████████[0m                        ]  51%[[32m██████████████████████████[0m                        ]  52%[[32m██████████████████████████[0m                        ]  53%[[32m███████████████████████████[0m                       ]  54%[[32m████████████████████████████[0m                      ]  55%[[32m████████████████████████████[0m                      ]  56%[[32m████████████████████████████[0m                      ]  57%[[32m█████████████████████████████[0m                     ]  58%[[32m██████████████████████████████[0m                    ]  59%[[32m██████████████████████████████[0m                    ]  60%[[32m██████████████████████████████[0m                    ]  61%[[32m███████████████████████████████[0m                   ]  62%[[32m████████████████████████████████[0m                  ]  63%[[32m████████████████████████████████[0m                  ]  64%[[32m████████████████████████████████[0m                  ]  65%[[32m█████████████████████████████████[0m                 ]  66%[[32m██████████████████████████████████[0m                ]  67%[[32m██████████████████████████████████[0m                ]  68%[[32m██████████████████████████████████[0m                ]  69%[[32m███████████████████████████████████[0m               ]  70%[[32m████████████████████████████████████[0m              ]  71%[[32m████████████████████████████████████[0m              ]  72%[[32m████████████████████████████████████[0m              ]  73%[[32m█████████████████████████████████████[0m             ]  74%[[32m██████████████████████████████████████[0m            ]  75%[[32m██████████████████████████████████████[0m            ]  76%[[32m██████████████████████████████████████[0m            ]  77%[[32m███████████████████████████████████████[0m           ]  78%[[32m████████████████████████████████████████[0m          ]  79%[[32m████████████████████████████████████████[0m          ]  80%[[32m████████████████████████████████████████[0m          ]  81%[[32m█████████████████████████████████████████[0m         ]  82%[[32m██████████████████████████████████████████[0m        ]  83%[[32m██████████████████████████████████████████[0m        ]  84%[[32m██████████████████████████████████████████[0m        ]  85%[[32m███████████████████████████████████████████[0m       ]  86%[[32m████████████████████████████████████████████[0m      ]  87%[[32m████████████████████████████████████████████[0m      ]  88%[[32m████████████████████████████████████████████[0m      ]  89%[[32m█████████████████████████████████████████████[0m     ]  90%[[32m██████████████████████████████████████████████[0m    ]  91%[[32m██████████████████████████████████████████████[0m    ]  92%[[32m██████████████████████████████████████████████[0m    ]  93%[[32m███████████████████████████████████████████████[0m   ]  94%[[32m████████████████████████████████████████████████[0m  ]  95%[[32m████████████████████████████████████████████████[0m  ]  96%[[32m████████████████████████████████████████████████[0m  ]  97%[[32m█████████████████████████████████████████████████[0m ]  98%[[32m██████████████████████████████████████████████████[0m]  99%[[32m██████████████████████████████████████████████████[0m] 100%                                                            

``` r
plot(jk.obj, xlab = "Variable Importance*100")
```

![](analysis0601_files/figure-markdown_github/unnamed-chunk-28-1.png)

## 6.8 Plot the estimated survival function for a simulated average patient

I use the random forest predictor created among all patients to estimate
survival of 3 new simulated patients. Their numeric features were
extracted from the medians in the population and categorical features
were extracted using the most frequent value in the population. The 3
patients are identical except for their treatment assignments. Their
estimated survival probability were plotted over time below. The plot
shows that the survival estimate for the 3 patients are almost
identical.

``` r
#create 3 patients data 
num_cols<-select_if(rf_1$xvar, is.numeric)
fac_cols<- select_if(rf_1$xvar, is.factor)
new_num<- data.frame(lapply(1:ncol(num_cols), function(i){
  median(num_cols[,i])
  }))
colnames(new_num)<- colnames(num_cols)
new_fac<- data.frame(lapply(1:ncol(fac_cols), function(i){
  which.max(table(fac_cols[,i]))
    }))
colnames(new_fac)<- colnames(fac_cols)
new<- cbind(new_num, new_fac)
new1<- new2<- new3<- new

#the three individuals only vary in terms of treatment 
new1[,which(rf_1$xvar.names=="treatment")]<- "hor_rad"
new2[,which(rf_1$xvar.names=="treatment")]<- "hor"
new3[,which(rf_1$xvar.names=="treatment")]<- "rad"
new_df<- rbind(new1, new2, new3)
y.pred<- predict(rf_1, newdata = new_df)
```

``` r
#plot survival estimates for the 3 patients 
plot_df<- data.frame(time= as.vector(y.pred$time.interest),
           new1=as.vector(y.pred$survival[1,]),
           new2=as.vector(y.pred$survival[2,]), 
           new3=as.vector(y.pred$survival[1,]))
plot_df %>% 
  ggplot(aes(x=time))+
  geom_line(aes(y=new1,colour= "hormone+radiation"), alpha=0.6,size=1.2)+
  geom_line(aes(y=new2, colour="hormone"), alpha=0.6, linetype="dotted", size=1.2)+ 
  geom_line(aes(y=new3, colour="radiation"),alpha=0.6, linetype="dashed", size=1.2)+
  xlab("Time (months) since study entry")+ 
  ylab("Survival probability")+
  ggtitle("Survival estimation for an average patient")+
  scale_x_continuous(breaks = c(12, 36, 60, 84, 108, 132, 180, 240, 300))+
  scale_color_manual(name="treatment", 
                     breaks=c("hormone+radiation", "hormone", "radiation"), 
                     values = c("red", "yellow", "green"))+
  theme_bw()
```

    ## Warning: Using `size` aesthetic for lines was deprecated in ggplot2 3.4.0.
    ## ℹ Please use `linewidth` instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

![](analysis0601_files/figure-markdown_github/unnamed-chunk-30-1.png)

# 7. Conclusion and recommendation

To summarize, I created predictors for survival among breast cancer patients who received different combos of hormone and radiation therapy. Because of the high dimension of the data that involves hundreds of gene expression features, I prioritize the predictive ability over interpretability (e.g., risk factors). The performance of these predictors are averagely in the range of 0.65-0.7, which is not excellent. However, the regularized cox hazard proportional models and the random survival forest all outperformed the traditional Cox proportional model by a lot, indicating of the potential of more accurate predictors with more observations. A major limitation of the study is the small number of observations relative to the number of features. Another limitation is that the time-to-event measurement do not start from a coherent time point (e.g., time of diagnosis of breast cancer), which makes it harder to interpret the results in terms of how to apply it to practice.

Due to limited data, further development is needed to enhance the accuracy and reliability of this model. One key step is ensuring that new clinical, mRNA, and mutation data are collected in a standardized manner, aligning with the current data collection plan. As more patient data become available, expanding the model to include individuals receiving other treatments, such as chemotherapy, would be beneficial. Additionally, continuous monitoring of key feature patterns, missing values, and predictive performance is essential. In particular, tumor stage has a high proportion of missing values, warranting further investigation into the reasons behind this gap and potential solutions for improvement. Given that the model aims to support clinical decision-making, the primary audience includes healthcare practitioners and potentially patients. Therefore, visualizing individual survival probability estimates over time, as demonstrated in this analysis, could be a valuable tool for interpretation and decision support.
