# creditmodel

`creditmodel` is a free and open source automated modeling R package designed to help model developers improve model development efficiency and enable many people with no background in data science to complete the modeling work in a short time.Let them focus more on the problem itself and allocate more time to decision-making.

`creditmodel` covers various tools such as data preprocessing, variable processing/derivation, variable screening/dimensionality reduction, modeling, data analysis, data visualization, model evaluation, strategy analysis, etc. It is a set of customized "core" tool kit for model developers.

`creditmodel` is suitable for machine learning automated modeling of classification targets, and is more suitable for the risk and marketing data of financial credit, e-commerce, and insurance with relatively high noise and low information content.

# Installation
```
# install.packages("creditmodel")
```
# Example
```
require(creditmodel)
if (!dir.exists("c:/test_model")) dir.create("c:/test_model")
setwd("c:/test_model")
#set parameters
LR.params = lr_params(
bins_control = list(bins_num = 8,bins_pct = 0.05, b_chi = 0.02, 
b_odds = 0.1,b_psi = 0.02,b_gb = 0.15,mono = 0.3,gb_psi = 0.05,kc = 1),
score_card = TRUE, cor_p = 0.7, iv_i = 0.02, psi_i = 0.1 )
XGB.params = xgb_params(nrounds = 10000, 
params = list(max.depth = 4, eta = 0.01, min_child_weight = 50, subsample = 0.5, colsample_bytree = 0.6, gamma = 0, max_delta_step = 1, eval_metric = "auc", objective = "binary:logistic"), early_stopping_rounds = 300)
#training model
Lending_model = training_model(
dat = lendingclub,
model_name = "lendingclub", target = "loan_status", occur_time = "issue_d",
ex_cols = c("last_credit_pull_d", "next_pymnt_d", "prncp|recoveries|rec_|funded_amnt|pymnt|fee$"),
obs_id = "id", prop = 0.7,
feature_filter = list(filter = c("IV", "PSI", "COR", "XGB"), cv_folds = 1, iv_cp = 0.02,
psi_cp = 0.1, cor_cp = 0.8,xgb_cp = 0, hopper = TRUE), algorithm = list("LR", "XGB"),
LR.params = LR.params, XGB.params = XGB.params,
parallel = FALSE,
save_pmml = FALSE,
 plot_show = FALSE,
seed = 46)
```
