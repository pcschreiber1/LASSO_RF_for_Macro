# Comp_Stat_Project
Project for the course of Computation statistics, summer semester 2021 MSc Economics University of Bonn

Content
 - Variable selection with Lasso
  - (Not really like bottom-up adding of variables to get optimal predictive model)
 - Variable selection with Random Forest (RF)
  - Why random forest and not boosting (who's penalty score looks similar to Lasso)
  - Why comparing linear and non-linear method?
    - common questions of interest like behaviour in collinearity?
      - both seem to **"bet" on sparsity**? [not entirely sure for RF]
    - or can they be competing methods? If so, for what type of applications
      - maybe also interesting to simply compare - since the latter one isn't used much in economics?
 - Simulation Study
    - Measures: MSE, comparison to Oracle rate (OLS with true DGP)
      - variable selection: retention frequency (Epprecht, 2013) [Genuer only use prediciton error]
 - Application
   
Covariance matrix $\rho^{i-j}$ or additional term for correlation of only some variables?

Genuer (2015): "RF is typically built using ntree = 2000 trees"

Is RF model selection also a shrinkage/regularization method or also a subset selection method?

Maybe not that important that they are the same method if they have other properties that connects the two: 
 - RF is supposed to have low variance, no? So are shrinkage methods!
    - So maybe the scenario where Lasso outperforms subset selection methods is also useful for RF? (Hastie et al 2020: Lasso outperforms forward step-wise selection if low signal-to-noise ratio (SNR))
 - Variable Importance scores are biased towards correlated variables (Genuer, 2015 - citing Strobl)
   - but not that helpful connection if both simply suffer from the same shortcoming (but maybe good for a side-note comparison)
 - Are LASSO and RF both embedded feature selection methods?
