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
#check missing after imputation (ranked)
missing<- data.frame(sapply(df_select_imp, function(x) sum(is.na(x)))) %>% 
  arrange(desc(sapply.df_select_imp..function.x..sum.is.na.x...))
head(missing, 10)


## 5.4 Feature selection based on correlation 

I dropped one of a pair of features if the pair correlation is over 0.7. I ended up with 497 features now. 

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

## 5.5 Stratify based on treatment 

There three main treatment types in this data, which are chemo, radiation, and hormone. Many patients have more than 1 treatment. Corresponding to my stated purpose, I decide to limit my data to patients who receive hormone and/or radiation therapy. I ended up with **1219 patients and 495** features. These patients are further stratified into hormone only, radiation only, and both for the purpose of creating different predictors.  

```{r}
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
table(df_select_imp$treatment)
#drop "chemotherapy", "hormone_therapy", "radio_therapy")) 
```

# 6. Model training and evaluation 

I generally followed the following steps: 
  
  - Explore survival curves of patients who had different treatment 
- Fit different models using the afromentioned algorithms and parameters tuning for patients of each of the 3 types of treatment 
- Evaluate the 10-fold cross-validated C-index
- Repeat the last two steps using the full data set of all patients who have had hormone and/or radiation therapy
- Compare performance

## 6.1 Survival curves of patients who received different treatment 

We can see that the unadjusted survival curves are significantly different among patients who receive different treatments. Echoing the concern about reduced efficacy of using both hormone and radiation therapy, we do see the worst survival for patients who receive both treatment. But this is likely confounded by other factors such as tumor stage and type.  

```{r}
cox1<- survfit(Surv(time, status)~ treatment, data=df_select_imp)
ggsurvplot(cox1, data= df_select_imp,
           title= "Survival curves of breast cancer patients by treatment group", 
           legend.title="",
           xlab= "Months since study entry",
           legend.labs=c("Hormone + Radiation Therapy", "Only Hormone Therapy", "Only Radation Therapy"),
           pval=T, conf.int = T, risk.table = T, tables.height = 0.2, tables.theme = theme_cleantable(), ggtheme=theme_bw(),
)
```

```{r}
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

```{r, warning=FALSE}
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

```{r, warning=FALSE}
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


```{r, warning=FALSE}
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

```{r, warning=FALSE}
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


```{r, warning=FALSE}
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


```{r, warning=FALSE}
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

```{r, warning=FALSE}
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


```{r, warning=FALSE}
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

```{r, warning=FALSE}
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


```{r, warning=FALSE}
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


```{r, warning=FALSE}
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

Based on the results, ridge predictor had the best performance in the predictor for patients who had hormone and radiation therapy. Lasso predictor had the best performance in the predictors for patients who only had hormone therapy. Random survival forest had the best performance for the predictor for patients who only had radiation therapy and overall patients. 

The overall performance is not that great, which might related to the fact that there are too many features vs. patients number. The performance using regularized regression and random survival forest are also close to each other. They all outperform the traditional cox proportional regression. 

```{r}
#patients had both hormone and radiation therapy
summary(eval1)
#patients only had hormone therapy
summary(eval2)
#patients only had radiation therapy
summary(eval3)
#among all patients
summary(eval_all)
```

# 6.7 Feature importance plot of the random survival forest predictor among all patients 

```{r}
set.seed(500)
rf_1<- rfsrc(Surv(time, status)~., data= df_select_imp, ntree=500, 
             samptype = "swr", importance=TRUE)
jk.obj<- subsample(rf_1)#default B=100 #of bootstrap
plot(jk.obj, xlab = "Variable Importance*100")
```

## 6.8 Plot the estimated survival function for a simulated average patient

I use the random forest predictor created among all patients to estimate survival of 3 new simulated patients. Their numeric features were extracted from the medians in the population and categorical features were extracted using the most frequent value in the population. The 3 patients are identical except for their treatment assignments. Their estimated survival probability were plotted over time below. The plot shows that the survival estimate for the 3 patients are almost identical. 

```{r}
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

```{r}
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

# 7. Conclusion 

To summarize, I created predictors for survival among breast cancer patients who received different combos of hormone and radiation therapy. Because of the high dimension of the data that involves hundreds of gene expression features, I prioritize the predictive ability over interpretability (e.g., risk factors). The performance of these predictors are averagely in the range of 0.65-0.7, which is not excellent. However, the regularized cox hazard proportional models and the random survival forest all outperformed the traditional cox proportional model by a lot, indicating of the potential of more accurate predictors with more observations.  

A major limitation of the study is the small number of observation relative to the number of features. Another limitation is that the time-to-event measurement seems do not start from a coherent time point (e.g.,time of diagnosis of breast cancer), which makes it harder to interpret the results in terms of how to apply it to practice.  

Because of limited data, this initial development needs further **implementation** steps to improve accuracy and reliability. One of the steps is to make sure that new clinical, mRNA and mutation data are measured in a standard way compared to the current data collection plan. As new patient data increase, it is ideal to expand this to patients who receive other types of treatment and combos (e.g., chemotherapy). We also need to constantly monitor the patterns of key features, missing values, and predictive performance over time. Specifically, tumor stage is a feature that has a lot of missing values. Further steps should takend to examine why this information is missing and how to improve. 

For the purpose of the predictors are to assist clinical decision making, I expect the audience to be practitioners in clinical settings or patients themselves. So, I think visualizing patients estimation of survival probability over time like I did above would be useful. I also think the eventual predictor should be able to automatically receive lab data (e.g., mRNA and mutation) for patient. It should also be an interactive tool where people can input and change fields if needed. However, it needs to be experimented whether imperfect predictions of one's lifespan would be beneficial or harmful overall.  
