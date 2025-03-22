# Predicting Survival of Breast Cancer Patients Using Gene Expression Data

# 1. Load libraries and data
rm(list = ls())
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

load("./processed_data/df_select_imp.rda")

# Check missing after imputation
missing <- data.frame(sapply(df_select_imp, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df_select_imp..function.x..sum.is.na.x...))
head(missing, 10)

# 2. Feature selection based on correlation
num_cols <- select_if(df_select_imp, is.numeric)
cor <- melt(round(cor(num_cols), 2)) %>% 
  arrange(desc(abs(value))) %>% 
  filter(Var1 != Var2)
big_cor <- cor %>% filter(abs(value) > 0.7)
drop_r <- as.vector(big_cor[, 1]) %>% unique()
df_select_imp <- df_select_imp[, !colnames(df_select_imp) %in% drop_r]
dim(df_select_imp)

# 3. Stratify based on treatment
df_select_imp <- df_select_imp %>% 
  mutate(treatment = case_when(
    chemotherapy == 1 & hormone_therapy == 1 & radio_therapy == 1 ~ "all",
    chemotherapy == 1 & hormone_therapy == 1 & radio_therapy == 0 ~ "che_hor",
    chemotherapy == 1 & hormone_therapy == 0 & radio_therapy == 1 ~ "che_rad", 
    chemotherapy == 1 & hormone_therapy == 0 & radio_therapy == 0 ~ "che", 
    chemotherapy == 0 & hormone_therapy == 1 & radio_therapy == 1 ~ "hor_rad", 
    chemotherapy == 0 & hormone_therapy == 0 & radio_therapy == 1 ~ "rad", 
    chemotherapy == 0 & hormone_therapy == 1 & radio_therapy == 0 ~ "hor", 
    chemotherapy == 0 & hormone_therapy == 0 & radio_therapy == 0 ~ "none"
  )) %>% 
  filter(treatment %in% c("hor", "hor_rad", "rad")) %>% 
  mutate(treatment = factor(treatment, levels = c("hor_rad", "hor", "rad"))) %>% 
  select(-chemotherapy, -hormone_therapy, -radio_therapy) %>% 
  rename(time = overall_survival_months, status = overall_survival)

dim(df_select_imp)
table(df_select_imp$treatment)

# 4. Survival curves by treatment group
cox1 <- survfit(Surv(time, status) ~ treatment, data = df_select_imp)
ggsurvplot(cox1, data = df_select_imp,
           title = "Survival curves of breast cancer patients by treatment group", 
           legend.title = "",
           xlab = "Months since study entry",
           legend.labs = c("Hormone + Radiation Therapy", "Only Hormone Therapy", "Only Radiation Therapy"),
           pval = TRUE, conf.int = TRUE, risk.table = TRUE,
           tables.height = 0.2, tables.theme = theme_cleantable(), ggtheme = theme_bw())

# 5. Split data by treatment groups
hor_rad <- df_select_imp %>% filter(treatment == "hor_rad") %>% select(-treatment)
hor <- df_select_imp %>% filter(treatment == "hor") %>% select(-treatment)
rad <- df_select_imp %>% filter(treatment == "rad") %>% select(-treatment)

# 6. Model evaluation helper function
run_models <- function(data, seed) {
  set.seed(seed)
  N <- nrow(data)
  V <- 10
  folds <- split(1:N, rep(1:V, length = N))
  eval <- matrix(NA, nrow = V, ncol = 5)
  colnames(eval) <- c("coxph", "lasso", "ridge", "el", "rf")
  
  Y <- data.matrix(data[c("time", "status")])
  X <- data.matrix(data[, !colnames(data) %in% c("time", "status")])
  
  for (v in 1:V) {
    train <- data[-folds[[v]], ]
    test <- data[folds[[v]], ]
    
    # Cox
    fit_cph <- coxph(Surv(time, status) ~ ., data = train)
    pred <- predict(fit_cph, newdata = test, type = "lp")
    eval[v, 1] <- concordance(Surv(time, status) ~ pred, data = test, reverse = TRUE)$concordance
    
    # Lasso
    fit_lasso <- cv.glmnet(x = X[-folds[[v]], ], y = Y[-folds[[v]], ], family = "cox", alpha = 1)
    pred <- predict(fit_lasso$glmnet.fit, newx = X[folds[[v]], ], s = fit_lasso$lambda.min, type = "link")
    eval[v, 2] <- concordance(Surv(time, status) ~ pred, data = data.frame(Y[folds[[v]], ], X[folds[[v]], ]), reverse = TRUE)$concordance
    
    # Ridge
    fit_ridge <- cv.glmnet(x = X[-folds[[v]], ], y = Y[-folds[[v]], ], family = "cox", alpha = 0)
    pred <- predict(fit_ridge$glmnet.fit, newx = X[folds[[v]], ], s = fit_ridge$lambda.min, type = "link")
    eval[v, 3] <- concordance(Surv(time, status) ~ pred, data = data.frame(Y[folds[[v]], ], X[folds[[v]], ]), reverse = TRUE)$concordance
    
    # Elastic Net
    fit_el <- cv.glmnet(x = X[-folds[[v]], ], y = Y[-folds[[v]], ], family = "cox", alpha = 0.5)
    pred <- predict(fit_el$glmnet.fit, newx = X[folds[[v]], ], s = fit_el$lambda.min, type = "link")
    eval[v, 4] <- concordance(Surv(time, status) ~ pred, data = data.frame(Y[folds[[v]], ], X[folds[[v]], ]), reverse = TRUE)$concordance
    
    # Random Forest
    fit_rf <- rfsrc(Surv(time, status) ~ ., data = train, ntree = 500, samptype = "swr")
    eval[v, 5] <- 1 - tail(predict(fit_rf, newdata = test)$err.rate, 1)
  }
  
  return(eval)
}

# 7. Run models by treatment group
eval1 <- run_models(hor_rad, seed = 777)
eval2 <- run_models(hor, seed = 555)
eval3 <- run_models(rad, seed = 333)
eval_all <- run_models(df_select_imp, seed = 111)

# 8. Model performance summary
summary(eval1)  # Hormone + Radiation
summary(eval2)  # Hormone only
summary(eval3)  # Radiation only
summary(eval_all)  # All patients

# 9. Random forest feature importance
set.seed(500)
rf_1 <- rfsrc(Surv(time, status) ~ ., data = df_select_imp, ntree = 500, samptype = "swr", importance = TRUE)
jk.obj <- subsample(rf_1)
plot(jk.obj, xlab = "Variable Importance*100")

# 10. Simulate new average patients and predict survival
num_cols <- select_if(rf_1$xvar, is.numeric)
fac_cols <- select_if(rf_1$xvar, is.factor)

new_num <- data.frame(lapply(num_cols, median))
new_fac <- data.frame(lapply(fac_cols, function(x) names(sort(table(x), decreasing = TRUE))[1]))
colnames(new_num) <- colnames(num_cols)
colnames(new_fac) <- colnames(fac_cols)

new <- cbind(new_num, new_fac)
new1 <- new2 <- new3 <- new
new1$treatment <- "hor_rad"
new2$treatment <- "hor"
new3$treatment <- "rad"

new_df <- rbind(new1, new2, new3)
y.pred <- predict(rf_1, newdata = new_df)

# 11. Plot survival estimates
plot_df <- data.frame(
  time = y.pred$time.interest,
  new1 = y.pred$survival[1, ],
  new2 = y.pred$survival[2, ],
  new3 = y.pred$survival[3, ]
)

plot_df %>% 
  ggplot(aes(x = time)) +
  geom_line(aes(y = new1, colour = "hormone+radiation"), alpha = 0.6, size = 1.2) +
  geom_line(aes(y = new2, colour = "hormone"), alpha = 0.6, linetype = "dotted", size = 1.2) + 
  geom_line(aes(y = new3, colour = "radiation"), alpha = 0.6, linetype = "dashed", size = 1.2) +
  xlab("Time (months) since study entry") + 
  ylab("Survival probability") +
  ggtitle("Survival estimation for an average patient") +
  scale_color_manual(name = "treatment", 
                     breaks = c("hormone+radiation", "hormone", "radiation"), 
                     values = c("red", "yellow", "green")) +
  theme_bw()
