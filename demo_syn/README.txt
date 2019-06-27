1. build traveltime table
cd tt_db
python taup_tt.py

2. synthetic test using a layer model
cd layermodel
python create_pk.py #create simulated P and S picks in './pk'
perl runREAL.pl #assocaite and locate these events using REAL and picks

3. synthetic test using homogenous/average velocities
cd homomodel
perl runREAL.pl #assocaite and locate these events using REAL and picks 
#picks are created based on a layer model, 
#but homogenous Vp and Vs are applied when using REAL
#traveltime table and searching settings are not necessary in this case

4. location comparision
cd plot
perl plotcatalogstation_layer.pl (GMT)
perl plotcatalogstation_homo.pl (GMT)

5. resolution analysis (example for synthetic test)
cd resolution
python create_pk.py (create P and S picks for one event)
perl runREAL.pl #resolution analysis only
plot_num.m (matlab) 

REAL usage:
see REAL warning (type REAL)
