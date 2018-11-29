This directory contains Matlab & XPP program files to implement the models and methods of sinoatrial node pacemaker cells coupled in a large network on 2D hexagonal grid.

Citation: Ly & Weinberg (2018). Analysis of Heterogeneous Cardiac Pacemaker Cell Models and Traveling Wave Dynamics. Journal of Theoretical Biology 459: 18-35.
DOI: 10.1016/j.jtbi.2018.09.023

We do not describe every m-file or sub-function because there are too many, rather we only describe the most pertinent files to reproduce our models and figures in the paper.

Sub-directories:    
runs_Hex/ — contains example scripts that were run for the full coupled Severi model of 271 cells; for 20s of simulation, takes a very long time (shorten time [line 58, tspan] to shorten run time)
!!must be moved to the main directory to run (relies on many .m and .mat files)!!
Naming convention follows: Hex[_Het/HomC/HomP]_[io/oi]_g[2/22/3/33/4]_s[1/2/3].m; Het is with all cell types (heterogeneous), HomC 
is a homogeneous network with all center cells, HomP is a homogeneous network with all peripheral cells; 
io=initial condition of network set to favor central wave generation (center->peripheral)
s1=gradient gap strength with peripheral having 15X larger strength than center.
s2=gradient gap strength with center having 15X larger strength than peripheral.
s3=constant gap strength throughout.
g2=0.25nS ([0.5/15,0.5]nS range of g_gap when gradient gap strength s1&s2), g22=0.625nS ([1.25/15,1.25]nS range with gradient gap), g3=1nS ([2/15,2]nS range with gradient gap), g33=1.5nS ([3/15,3]nS range with gradient gap), g4=2nS ([4/15,4]nS range with gradient gap)

thry/ - theory; coupled oscillator models, contains files to implement the phase reduction model.  See 1README.txt in that directory


Main Files for Severi Model:     
SA_fcn_Het.m — the main ODE file (right-hand-side) to implement the Severi Model, with 31 dynamic variables.  
Input parameters: (t,y) usual ODE, t=time, y=state variable (31*Nz sized vector). Nz=total # pacemaker cells, G_cup=gap junctional coupling matrix; 
the next 12 inputs are intrinsic cellular parameters to give heterogeneity, each is an Nz x 1 vector: 
C,Lcell,Rcell,gfNa,gfK,PCaL,PCaT,gKr,gKs,gNa,INaKmax,KNaCa
SA_fun_Interm.m — for running a single Severi cell model, much like Severi et al. ’12 BUT allows heterog in the 12 parameters we varied

Mat Files:     
disi.mat — contains variable isi_L that is an 18 x 1 vector of the cycle lengths (CL) or interspike internal (ISI), in seconds.  Calculated by running the model and using the function get_isiL.m AFTER running SA_Model_Int_Plots.m (SA_fun_Interm.m) and inputting (tm,Volt) and threshold for spiking
Parms_c2p.mat — contains heterogeneous parameter values of all 18 cell types; variables: C,Lcell,Rcell,gfNa,gfK,PCaL,PCaT,gKr,gKs,gNa,INaKmax,KNaCa, 
each of size 18 x 1.  Also has scl_vl and iv_scl which are 18 x 1 vectors to denote how parameters were linearly interpolated between Center 
and Peripheral (except INaKmax, KNaCa, which were manually altered to insure a stable limit cycle).
LimC_IC_all.mat — contains 18 x 1 cell variable, where each cell IC_cl{j} is a matrix with 31 columns denoting the values of the 31 variables 
over 1 cycle (uncoupled) of that particular (j) cell type (spike to spike).  This is used for setting initial conditions, to set how near or far cell is 
from spiking; the first row (1 x 31) has state variable values furthest away from spiking, last row (1 x 31) has state variable values closest to spiking.
* NOTE: although there are 18 possible cell types, we only use 9, index via: 1,3,5,7,9,11,13,15,18

Other m-file functions:
createHexGrid.m — script that returns the spatial location (X_l,Y_l) of all Nz cells on hexagonal grid; after running saves variables (Nz,X_l,Y_l,ind_j) in  H_dataHexGrd.mat [see function comments for further details].
nn_coupVryStrHex.m — function to generate the coupling matrix G_cup.  Inputs: (X_l,Y_l) spatial location, gp_lim=Euclidean distance for cells to be coupled 
(too small=>uncoupled, too large=>all-to-all coupling), strCent (gradient C->P?), ind_j (identify cell types); in this paper we set so ONLY nearest neighbor coupling, with gp_lim=1.05.
svPhaseIC.m - gets the initial distribution (IO or OI) and transforms to phase variables [0,1] for all 3 regimes (Het, HomC, Home), saves result in ~/thry/Phase_IC.mat.

show_cellsConnectHex.m - a script to show network configuration/connect (relies on  H_dataHexGrd.mat)

Generate Main Figures:

Figure 1
A) Run SA_Model_Int_Plots.m, which REQUIRES user input.  Type in: 1,3,5,7,9,11,13,15,18 to get the 9 cell types used in the paper.  After running 9 times, have voltage trajectories lined up with all 9 different cell types in Figure 39.  Zoom-in around 1 spike (time 0).  Script calls on SA_fun_Interm.m.
B) Same data as in A), just displayed differently.
C) Data saved in disi.mat (generated from get_isiL.m, see above), in variable isi_L which is an 18 x 1 vector.  Again, we only use entries: 1,3,5,7,9,11,13,15,18 .  Inset (frequencies) is 1/isi_L.

Figure 2
A) Run script plotVolt_traj.m and enter inputs in sequence: 2, 2, 3.
B) Run script plotVolt_traj.m and enter inputs in sequence: 3, 2, 3.
D) Run script plotVolt_traj.m and enter inputs in sequence: 1, 2, 3.
C) Left: run script show_cellsConnectHex.m; Right: show_cellsInitC_Hex.m

Figure 3
A) Same data as in Fig 1 C (disi.mat) but scaled by isi_L(1); Freq=isi_L(1)/(isi_L(j)), for j=1,3,5,7,9,11,13,15,18
B) In ~/thry/, run plot_P_Hs.m
C) In ~/thry/, run plot_circ_chains.m

Figure 4
Use the script ~/thry/plotHfcn_suffCondit_hex.m and input the corresponding network type (Het, Homog C, Homog P), 
identifier for coupling strength, and coupling type.  Relies on mat files: sHex_H[1/2/3]_g[2/22/3/33/4]_s[1/2/3].mat, created from phaseLags_h_MatchHex.m (see description of that file in ~/thry/1README.txt)
See Fig S3 for similar plots, different parameters.

Figure 5
Close all figures in Matlab first; 
For the solid lines (theoretical freq from reduced phase model), run ~/thry/plot_Freqcs.m; again, relies on the same mat files as in Figure 4, sHex_H[1/2/3]_g[2/22/3/33/4]_s[1/2/3].mat
For the stars (long simulations), run plots_Freq.m in main directory (creates mat file Freq_fullModel.mat, which is used later).  !!ASSUMING solid lines from theory are in Figures 1, 2, 3!!

Figure 6
Panels A, C, E has 2 sets of curves; run ~/thry/plot_TimeToSS_eig.m FIRST, keeping 3 figures open, THEN run plots_timeSS.m (!!! assumes that 3 figures are in Figure 1, Figure 2, Figure 3!!!)
Panels B, D, F (left column) are plotted with ~/thry/plot_CharTimeScale.m; again, relies on the same mat files as in Figure 4, sHex_H[1/2/3]_g{2/22/3/33/4]_s[1/2/3].mat  
1) The time to steady state from the phase model is plotted via ~/thry/plot_TimeToSS_eig.m, which relies on same sHex_H[].mat functions; plots gray curves that use the entire eigenfunction expansion (rather than just the largest negative real big).  Running this will create a mat file that saves times in  d_tmeSS_IO.mat [used in Figures 7, 8, 9]. 
2) The time steady state from the full model is plotted via plots_timeSS.m, which loads simulation results in mat files (dhex_het_[], dhexhomC_[], dhexhomP[]), Freq_fullModel.mat (from plots_Freq.m), and ~/thry/svGavgCompHex.mat — after running this, will create file dTimeSS_all.mat to save results [used in Figures 7, 8, 9]

Figures 7-9
All generated from plotVolt_wTransTime.m, which loads simulation results in mat files (dhex_het_[], dhexhomC_[], dhexhomP[]), and d_tmeSS_IO.mat (transient times from phase model) and dTimeSS_all (transient times from the full large-scale models).  Requires user input to show particular figure.

Figure 10
Obtained from a very LONG run of data similar to in Figs 4,5,6 [combining sim_phaseOscill_MatchHex.m, get_y0_fromSims.m, phaseLags_h_MatchHex.m into long automated program].  

Figure 11
A), B), and C) generated by the script: show_nonTW.m
D), E), and F) generated by 

Figure S1 
See how Figure 4 was generated.

Figure S2
Run script ~/thry/plots_timeSS.m

Figure S3-S5
See how Figures 7-9 generated.

Figure S6-S9
Same as prior figures but with larger and smaller g_gap; use scripts and results from Figure 10, re-run larger model with phase reduced model parameters
