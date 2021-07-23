# Rapid Earthquake Association and Location (REAL)
The REAL package contains following main files:
1) REAL: code for Rapid Earthquake Association and Location 
2) VELEST: code for 1-D inversion of velocities and hypocenter locations
3) domo_syn: demo for synthetic tests
4) demo_real: demo for real data (data downloading-> picking -> associaiton -> location -> relocation)

A more comprehensive workflow will be available at https://github.com/Dal-mzhang/LOC-FLOW

# Usage:
See user guide and README files
(raw data -> high-resolution earthquake locations)

# Introduction:
REAL (Rapid Earthquake Association and Location) associates arrivals of different seismic phases and locates seismic events primarily through counting the number of P and S picks and secondarily from traveltime residuals. A group of picks are associated with a particular earthquake if there are enough picks within the theoretical traveltime windows. The location is determined to be at the grid point with most picks. If multiple locations have the same maximum number of picks, the grid point with smallest traveltime residual is selected. We refine seismic locations using a least-squares location method (VELEST) and a high-precision relative location method (hypoDD).

# References:
1) Zhang, M., W.L. Ellsworth, and G.C. Beroza. Rapid Earthquake Association and Location, Seismol. Res. Lett., 90.6, 2276-2284, 2019, https://doi.org/10.1785/0220190052
2) Liu M., Zhang M., Zhu W., Ellsworth W. and Li H. Rapid Characterization of the July 2019 Ridgecrest, California Earthquake Sequence from Raw Seismic Data using Machine Learning Phase Picker. Geophysical Research Letters, 47(4), e2019GL086189, 2020, https://doi.org/10.1029/2019GL086189
3) Wang R., Schmandt B., Zhang M., Glasgow M., Kiser E., Rysanek S. and Stairs R. Injection-induced earthquakes on complex fault zones of the Raton Basin illuminated by machine-learning phase picker and dense nodal array. Geophysical Research Letters, 47(14): e2020GL088168, 2020, https://doi.org/10.1029/2020GL088168
4) Kissling, E., W.L. Ellsworth, D. Eberhart-Phillips, and U. Kradolfer: Initial reference models in local earthquake tomography, J. Geophys. Res., 99, 19635-19646, 1994, https://doi.org/10.1029/93JB03138
5) Waldhauser F. and W.L. Ellsworth, A double-difference earthquake location algorithm: Method and application to the northern Hayward fault, Bull. Seism. Soc. Am., 90, 1353-1368, 2000, https://doi.org/10.1785/0120000006

# Author:
Miao Zhang, Dalhousie University, miao.zhang@dal.ca

# Versions:
REAL1.0, June  27, 2019

The codes are improvded over time, see changes in the Modified_History.
