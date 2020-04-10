# geostat

geostat aims to analyze on-farm research data.

'likfit2' fits the anisotropic model (sum-metric model) with exponential covariance functions.
'variog2' evaluate the experimental variograms. 

Two dimentions are need (e.g. direction of tractors' travel (dist1) and perpendicular to travel (dist2)). It is implemented by MultiStart function. Eight values [nugget sill1 sill2 sill3 rho1 rho2 rho3 alpha] should be optimized. Initial values (x0), the number of random runs (Nrun), lower and upper bounds for values (lower, upper) should be specified. Two estimators (ML and REML) are available. REML will be used if a vale for REML is 1.Ten runs may take 3min.
