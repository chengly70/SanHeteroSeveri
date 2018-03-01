
Describing sub-directory ~/thry/ that has model/methods for 
reduced 1D chain of oscillators with cut-end (nearest neighbor coupling) approximation to the full 2D coupled Severi Model


The H functions are calculated in XPP via the file: sa_Int_jth.ode; results are saved Hfs.mat

Primary M-files:

sim_phaseOscill_MatchHex.m — simulate the reduced phase models with attributes from full Severi network with proper H-fcn & 
frequencies, see if settles to traveling wave.
Lines 61-68 show how we systematically calculate the effective H functions from the full Severi network

get_y0_fromSims.m — getting the ‘correct’ y0 for phaseLags_h_MatchHex routine, AFTER running sim_phaseOscill_MatchHex.m; 
user MUST input where periodic steady-state starts (view figures from sim_phaseOscill_MatchHex.m)

phaseLags_h_MatchHex.m — main file to numerically calculate traveling wave solution, AFTER running sim_phaseOscill_MatchHex.m and get_y0_fromSims.m
saves results in ss_H[]_g[]_s[].mat

eqn_distr.m — the objective function for the periodic steady-state, using H_j, but allowing cell distributions via vit (diff # of cell types)
Used in phaseLags_h_MatchHex.m

plot_Freqs.m - plots the frequency from the phase model; relying 
on numerical solutions to periodic steady state eqns (data from Fig 4, see phaseLags_h_MatchHex.m)

plot_CharTimeScale.m — Gets the characteristic time scale (the eigenvalue with largest negative real part); again relies on data from Fig 4 to calculate the fundamental matrix of how perturbations evolve.

Other M-files, post-processing:

plot_P_Hs.m — script to plot the Hfunctions (& PRCs) from XPP
plotHfcn_suffCondit_hex.m — script to plot all Hfunctions to verify sufficient condition from Ermentrout Theorem.
User MUST specify which network, coupling, etc.
plot_Freqcs.m — Script to plot the frequency of the traveling wave (reduced phase model) with physiological units (1/s).
Relies on svGavgCompHex.mat, which calculates the average g_gap values for a given (full 2D) Severi network.  svGavgCompHex.mat is from running calc_gAvg_fromFullHex.m, 
which loads H_dataHexGrd.mat


Mat-Files:
Hfs.mat — the computed interaction (H) functions from XPP, for all 18 cells assuming coupled to self.
svGavgCompHex.mat - average g_gap from full 2D Severi; 3 rows representing s1, s2, s3 coupling; 5 columns with 
g2, g22, g3, g33, g4 (same convention, see below or other 1README.txt file)
PRCs.mat — numerically calculated PRCs for all different cell types

Results from running phaseLags_h_MatchHex.m:
sHex_H[1/2/3]_g[2/22/3/33/4]_s[1/2/3].mat, 
Naming convention for ss_H*.mat follows:
H1=Het is a ‘uniform’ distribution of heterogeneous cells. 
H2=HomC is a homogeneous network with all center cells. 
H3=HomP is a homogeneous network with all peripheral cells. 
s1=gradient gap strength with peripheral having 15X larger strength than center.
s2=gradient gap strength with center having 15X larger strength than peripheral.
s3=constant gap strength throughout. 
g2=0.25nS (0.5nS largest g_gap when gradient gap strength s1&s2). 
g22=0.625nS (1.25nS with gradient gap). 
g3=1nS (2nS with gradient gap). 
g33=1.5nS (3nS with gradient gap). 
g4=2nS (4nS with gradient gap).
These attributes refer to full large-scale Severi model, which we systematically map to the 1D chain of phase oscillators.
