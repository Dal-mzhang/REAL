# REAL
Rapid Earthquake Association and Location
                   
                   REAL1.0    2019/06/11
              https://github.com/Dal-mzhang/REAL
                Miao Zhang, Dalhousie University 
                    E-mail: miao.zhang@dal.ca

The package is free and distributed in the hope that it will be useful. 
You can use it and/or modify it under the terms of the GNU General Public License 
for NON-COMMERCIAL purposes.

It does not come with any warranties, nor is it guaranteed to work on your computer. 
The user assumes full responsibility for the use of this system. 
We are NOT responsible for any damage that may follow from correct or incorrect 
use of this package.

Permission is granted to make and distribute verbatim copies of this license provided 
that the copyright notice and these paragraphs are preserved on all copies.

1. Content of the Package
    The REAL package contains following main files:
    domo_syn: demo for synthetic tests
    deom_real: demo for real data
    REAL: code for Rapid Earthquake Association and Location 
    VELEST: code for 1-D inversion of velocities and hypocenter locations

2. Installation
    Please see Makefile in each directory.

3. Usage
    Please see README in each directory. type 'REAL'.

4. Workflow
    Please see Fig. 2 (Flow diagram of the REAL method) in the reference paper.

5. Reference
   Zhang, Ellsworth and Beroza, Rapid Earthquake Association and Location, SRL (in revision)

6. Release History
    REAL Version 1.0     2019/06/11

7. Requirement (optional) 
    VELEST (http://www.seg.ethz.ch/software/velest.html, for location refinement, already provided here)
    hypoDD (https://www.ldeo.columbia.edu/~felixw/hypoDD.html, for location refinement)
    ObsPy (https://github.com/obspy/obspy, for data downloading, processing, and picking only)
    SAC (http://ds.iris.edu/ds/nodes/dmc/software/downloads/sac, for data processing only)
    GMT (https://www.soest.hawaii.edu/gmt/, for figure plotting only)
    Matlab (https://www.mathworks.com/products/matlab.html, for figure plotting only)

8. Attentions
    (1) REAL only needs seismic picks. Picks can be obainted STA/LTA (provided) or other advanced methods (e.g., maching-learning)
    (2) REAL doesn't depend on any other softwares but requires VELEST and hypoDD to further refine earthquake locations
    (3) Parameters settings are empirical. You may change them for your case.
