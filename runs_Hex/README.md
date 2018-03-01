Sub-directory:
runs_Hex/ — contains example scripts that were run for the full coupled Severi model of 271 cells; for 20s of simulation, takes a very long time (shorten time [line 58, tspan] to shorten run time)
!!must be moved to the main directory to run (relies on many .m and .mat files)!!
Naming convention follows: Hex[_Het/HomC/HomP]_[io/oi]_g[2/22/3/33/4]_s[1/2/3].m; Het is with all cell types (heterogeneous), HomC 
is a homogeneous network with all center cells, HomP is a homogeneous network with all peripheral cells; 
io=initial condition of network set to favor central wave generation (center->peripheral)
s1=gradient gap strength with peripheral having 15X larger strength than center.
s2=gradient gap strength with center having 15X larger strength than peripheral.
s3=constant gap strength throughout.
g2=0.25nS ([0.5/15,0.5]nS range of g_gap when gradient gap strength s1&s2), g22=0.625nS ([1.25/15,1.25]nS range with gradient gap), g3=1nS ([2/15,2]nS range with gradient gap), g33=1.5nS ([3/15,3]nS range with gradient gap), g4=2nS ([4/15,4]nS range with gradient gap)
