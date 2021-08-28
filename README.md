# Comp_Stat_Project
--- 
Project for the course in Computational Statistics | Summer 2021, M.Sc. Economics, Bonn University | [Philipp Schreiber](https://github.com/pcschreiber1)

# Variable Selection with High-Dimensional Low-Quality Data
## A comparison of LASSO and Random Forest. <a class="tocSkip">   
    
    
---

In this notebook, I contribute to the analysis of variable selection under high-dimensional low quality data by comparing the performance of the LASSO and relaxed Lasso to another method designed for efficient varaince reduction: Random Forests (RF). RF-based variable selection procedures have become especially popular in bio-medical sciences, but, so-far, are still uncommon in many strands of economics. Here we apply the data-driven technique of Genuer et al. (2010), which is implemented in the `VSURF` package.  Three different simulation studies are conducted which closely emulate real-world macro-economic data and each illustrate the methodologies' behaviour under different challenges. Importantly, given the liability of direct comparisons of parametric and non-parametric methods to the underlying data generating process (DGP), we here focus on the techniques' relative performance under different levels of noise clarity. The discussion is complemented by an application to Sala-I-Martin's (1997b) famous "millions" data set, which has frequently been used to showcase variable selection strategies for macro-economic growth models.