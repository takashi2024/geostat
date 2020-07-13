## geostat v0.1.0
'geostat' aims to analyze on-farm research data. 'geostat' is developed using Matlab 2020a. The effects of treatment (e.g. fertilizer application) or environental variables (soil properties) on crop yield or quality can be evaluated using a spatial linear mixed model. The spatial linear mixed model can account for spatial autocorrelation in a response variable which leads to unreliable inferences. Yield monitor or remotely senesed data are assumed to be used as a response variable.

## Requirement
* Optimization Toolbox
* Global Optimization Toolbox
* Parallel Computing Toolbox (optional)

Environments under Matlab 2020a are tested.

## How to use
Please run each section code in 'demo.m' to learn how to use following functions. Three dataset ('Field1–3.csv') are availble for this demo.

* '__likfit__' fits the isotropic model.<br>
Exponential or Matern covariance functions are available. Three values [nugget sill rho] will be optimized for the exponential covariance function. Initial values (x0), the number of random runs (Nrun*), lower and upper bounds for values (lower, upper) should be specified. Two estimators (ML and REML) are available. REML will be used if a value for REML is 1. This function is basically originated from R package 'geoR'.

* '__variog__' evaluates omni-directional experimental variograms.

* '__likfit2__' fits the anisotropic model (sum-metric model).<br>
Two dimentions are needed (e.g. direction of tractors' travel (x) and perpendicular to travel (y)) for the parameter estimations. Sum-metric model has been widely used for space-time modelling in the previous literatures, and it is developed for anisotropic modelling.  Exponential or Matern covariance functions are available currently. Eight values [nugget sill1 sill2 sill3 rho1 rho2 rho3 alpha] should be optimized for the exponential convariance function. Initial values (x0), the number of random runs (Nrun*), lower and upper bounds for values (lower, upper) should be specified. Two estimators (ML and REML) are available. REML will be used if a vale for REML is 1.

* '__variog2__' evaluates two directional experimental variograms.

*To avoid the local minimum problem in minimizing negative log-likelihood, initial parameter sets will be randomly chosen and fitted for Nrun times by Global Optimization Toolbox.

## Code sharing policy
The code is available without limitations for any purposes. However, if a publication, conference paper, book, poster, or any other publications as an outcome of research work, reference is required. In case a help is needed in using the code, or adapting the code to a specific dataset, please do not hesitate to contact me. Takashi Tanaka, 11/04/2020

## Author
* Takashi S. T. Tanaka
* Faculty of Applied Biological Sciences, Gifu Univresity
* Twitter: https://twitter.com/SonamTashi331
* Researchgate: https://www.researchgate.net/profile/Takashi_Tanaka8

## Notice
The manuscript which uses this source code is in preparation. Preprint version is available from https://arxiv.org/abs/2004.12741.
