# geostat
geostat aims to analyze on-farm research data. Yield monitor or remotely senesed data are assumed to used as a response variable.


'likfit2' fits the anisotropic model (sum-metric model) with exponential covariance functions.
'variog2' evaluate the experimental variograms. 

# How to use
Two dimentions are need (e.g. direction of tractors' travel (x) and perpendicular to travel (y)). It is implemented by MultiStart function. Eight values [nugget sill1 sill2 sill3 rho1 rho2 rho3 alpha] should be optimized for exponential convariance model. Initial values (x0), the number of random runs (Nrun), lower and upper bounds for values (lower, upper) should be specified. Two estimators (ML and REML) are available. REML will be used if a vale for REML is 1.

# Code sharing policy
The code could be used without limitations for any purposes. However, if a publication, conference paper, book, poster, or any other publications as an outcome of such a work, reference is required. In case a help is needed in using the code, or adapting the code to a specific dataset, please do not hesitate to contact me. Takashi Tanaka, XX/04/2020
