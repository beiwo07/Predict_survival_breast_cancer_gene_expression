# Predicting Survival of Breast Cancer Patients Using High-dimensional Gene Expression Information

# Clear environment
rm(list = ls())

# Load libraries
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

# Read data
df <- read.csv("raw_data/METABRIC_RNA_Mutation 2.csv")
dim(df)

# Data summary
# Binary overall survival status
table(df$overall_survival)

# Age 
summary(df$age_at_diagnosis)

# Cancer type 
table(df$cancer_type_detailed)

# Treatment 
table(df$chemotherapy)
table(df$hormone_therapy)
table(df$radio_therapy)
table(df$chemotherapy, df$hormone_therapy, df$radio_therapy)

# Number of mutations 
table(df$mutation_count)

# Clinical variables overview
clinical <- df[1:31]
sapply(clinical, function(x) unique(x))

# Missing value summary
missing <- data.frame(sapply(df, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df..function.x..sum.is.na.x...))  
head(missing, 10)

# Check missing values in outcomes
table(is.na(df$overall_survival))
table(is.na(df$overall_survival_months))

# Replace empty strings with NA
df <- df %>% 
  replace_with_na_all(condition = ~.x == "")
clinical <- clinical %>% 
  replace_with_na_all(condition = ~.x == "")

# Drop non-useful features
drop <- c("cancer_type", "cohort", "patient_id", "er_status_measured_by_ihc", 
          "her2_status_measured_by_snp6", "nottingham_prognostic_index", 
          "death_from_cancer", "overall_survival", "overall_survival_months")
outcome <- df[, colnames(df) %in% c("overall_survival", "overall_survival_months")]
df_select <- df[, !colnames(df) %in% drop]

# Drop sparse mutation features
df_select <- df_select %>% select(-ends_with("_mut"))

# Check dimensions
dim(df_select)
dim(outcome)

# Check missingness again
missing <- data.frame(sapply(df_select, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df_select..function.x..sum.is.na.x...)) 
head(missing, 15)

# Convert character features to factor before imputation
df_select[sapply(df_select, is.character)] <- lapply(df_select[sapply(df_select, is.character)], 
                                                     as.factor)

# Multiple Imputation by Chained Equations (optional to evaluate)
# set.seed(777)
# imp <- mice(df_select[, !colnames(df_select) %in% c("overall_survival_months", "overall_survival")], 
#             m = 3, method = "pmm")
# df_select_imp <- complete(imp)

# Add outcome back
# df_select_imp <- data.frame(df_select_imp, outcome)

# Save processed data
# save(df_select_imp, file = "/processed_data/df_select_imp.rda")
