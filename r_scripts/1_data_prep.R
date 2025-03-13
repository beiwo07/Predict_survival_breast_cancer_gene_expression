library(tidyverse)
library(ggplot2)
library(naniar)
library(mice)
library(reshape2)
library(survival)
library(glmnet)
library(survminer)
library(randomForestSRC)
library(iml)
library(magrittr)

#original dataset dimension
df<- read.csv(".raw_data/METABRIC_RNA_Mutation 2.csv") 
dim(df)

# 4.Feature engenieering and exploratory data analysis

## 4.1 Data dimension & key variable distributions 

# The orginal dataset contains 1904 observations and 693 features. Among the features, there are 31 clinical variables. There are 331 mRNA level z-scores (numeric) for 331 genes and mutation (binary) for 175 genes. Data on mutation of genes are very sparse, so that I decided to just focus on clinical features and mRNA data. 

# My exploration indicates that 801 deaths occurred during the follow up period. The average age of the cohort if 61. Patients in this cohort receive three types of therapies and combinations of them. The data also shows that very few patients have mutation on more than 15 genes. I also found that **missing values** were only partially coded as NA. Some of them are empty cells. 

#binary overall survival status
table(df$overall_survival)
#age 
summary(df$age_at_diagnosis)
#cancer type 
table(df$cancer_type_detailed)
#treatment 
table(df$chemotherapy)
table(df$hormone_therapy)
table(df$radio_therapy)
table(df$chemotherapy, df$hormone_therapy, df$radio_therapy)
#number of mutation 
table(df$mutation_count)

#take a closer look at clinical variables 
clinical<- df[1:31]
sapply(clinical, function(x) unique(x))

#missing values ranked; not all are correctly represented
missing<- data.frame(sapply(df, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df..function.x..sum.is.na.x...))  
head(missing, 10)

# 5. Data pre-processing 

In my data pre-processing, I follow the steps below: 
  
  - Explore missing values and correctly code missing as NAs
- Manually drop features that will not be useful for this task 
- Multiple imputation by chained equation
- Create categorical variables to stratify patients based on treatment 

## 5.1 Explore missing values and replace missing with NAs

#there isn't missing values in the outcome 
table(is.na(df$overall_survival))
table(is.na(df$overall_survival_months))

#replace empty cells with NA 
df<- df %>% 
  replace_with_na_all(condition= ~.x=="")
#a subset of clinical variables for viewing purpose
clinical<- clinical %>% 
  replace_with_na_all(condition= ~.x=="")

## 5.2 Drop clinical features that will not be useful 

#clinical feature list
str(clinical)

#drop redundant, not useful, outcome features.. 
drop<-c("cancer_type", "cohort","patient_id", "er_status_measured_by_ihc", "her2_status_measured_by_snp6","nottingham_prognostic_index","death_from_cancer", "overall_survival", "overall_survival_months")
outcome<- df[, colnames(df) %in% c("overall_survival", "overall_survival_months")]
df_select<-df[, !colnames(df) %in% drop]

#drop mutations in genes because too sparse 
df_select<- df_select %>%
  select(-ends_with("_mut"))

#dims after dropping 
dim(df_select)
dim(outcome)

## 5.3 Multiple imupatation 

Now that all missing are correctly coded, I first rank the features by numbers of missingness. I can see that the missing values are all on the clinical features. I use multiple imputation with chained equations depending on the type of the features. The data table below looks like there is no missingness. 

The imputation process is muted here because I saved the imputed dataset to save time. But the code is shown here. 

#rank by missing #
missing<- data.frame(sapply(df_select, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df_select..function.x..sum.is.na.x...)) 
head(missing, 15)

#convert all character features to factor for imputation
df_select[sapply(df_select, is.character)] <- lapply(df_select[sapply(df_select, is.character)], 
                                                     as.factor)
#imputation (excluding the outcome features)
set.seed(777)
imp<- mice(df_select[, !colnames(df_select) %in% c("overall_survival_months","overall_survival")], m=3, method = "pmm")
df_select_imp<- complete(imp)

#add outcome back
df_select_imp<- data.frame(df_select_imp, outcome)

#save this imputed data so that I don't have to run again 
save(df_select_imp, file="./processed_data/df_select_imp.rda")
