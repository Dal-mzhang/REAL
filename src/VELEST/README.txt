original velest could be downloaded at http://www.seg.ethz.ch/software/velest.html
I met some issues when I complied it on MAC. Here I made some changes to fix them.
Main changes by M. Zhang:
1. changed some minor format issue so that it can be complied.
2. changed the station 4 char to  6 char
3. changed the output velocity format (isingle=0),now you can directly use it as input velocity
4. fixed a bug (in previous version,the last pick always missing) 

Note: you may change/increase ieq in vel_com.f in your case. The default value is 658 (number of events). 

UserGuide: https://www.ethz.ch/content/dam/ethz/special-interest/erdw/geophysics/seismology-and-geodynamics-dam/document/velest_guide.pdf
Reference: Kissling, E., W.L. Ellsworth, D. Eberhart-Phillips, and U. Kradolfer: Initial reference models in local earthquake tomography, J. Geophys. Res., 99, 19635-19646, 1994
