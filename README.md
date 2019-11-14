# Rapid Earthquake Association and Location (REAL)
The REAL package contains following main files:
1) REAL: code for Rapid Earthquake Association and Location 
2) VELEST: code for 1-D inversion of velocities and hypocenter locations
3) domo_syn: demo for synthetic tests
4) demo_real: demo for real data

# Usage:
See user guide and README files
(raw data -> high-resolution earthquake locations)

# Introduction:
REAL (Rapid Earthquake Association and Location) associates arrivals of different seismic phases and locates seismic events primarily through counting the number of P and S picks and secondarily from traveltime residuals. A group of picks are associated with a particular earthquake if there are enough picks within the theoretical traveltime windows. The location is determined to be at the grid point with most picks. If multiple locations have the same maximum number of picks, the grid point with smallest traveltime residual is selected. We refine seismic locations using a least-squares location method (VELEST) and a high-precision relative location method (hypoDD).

# References:
1) Zhang, M., W.L. Ellsworth, and G.C. Beroza. Rapid Earthquake Association and Location, Seismol. Res. Lett., 90.6, 2276-2284, 2019.
2) Kissling, E., W.L. Ellsworth, D. Eberhart-Phillips, and U. Kradolfer: Initial reference models in local earthquake tomography, J. Geophys. Res., 99, 19635-19646, 1994.
3) Waldhauser F. and W.L. Ellsworth, A double-difference earthquake location algorithm: Method and application to the northern Hayward fault, Bull. Seism. Soc. Am., 90, 1353-1368, 2000.

# Author:
Miao Zhang, Dalhousie University, miao.zhang@dal.ca

# Versions:
REAL1.0, June  27, 2019

REAL1.1, Sept. 13, 2019

REAL 1.2, Nov. 14, 2019
