tools to make plots are here

plots at "model" 
dsigma_dt and sigma made by cross_section.ipynb to 
flux made by flux.C

Harry Lee's ed differential crossection model is on a grid of Eg and t with limited range
interpolation and extrapolation are needed to cover full kinematic range
jpsi_t.C will output the kinematic ranges
Eg      t0              t1              t at 28deg of SoLID max theta acceptance
6       -1.86728        -6.68353        -3.03244
7       -0.952465       -11.4521        -3.49256
8       -0.620395       -15.6147        -4.24786
9       -0.444399       -19.605         -5.07979
10      -0.336626       -23.5154        -5.94409
12      -0.214367       -31.2184        -7.71494
14      -0.149221       -38.8427        -9.51003
16      -0.110103       -46.4269        -11.3152
18      -0.0846834      -53.988         -13.1251

plots made by plot.C at
"ed_dvmp_8.8GeV_photoproduction_extranear" where extrapolation is equal to nearest points, this overestimates at high t
"ed_dvmp_8.8GeV_photoproduction_extrano" where no extrapolation is done and dsigma_dt is 0 for the outside kinematic range

