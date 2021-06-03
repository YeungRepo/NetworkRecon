### Subrepo manifest: 

*.mat files contain experimental data from the in vitro repressilator validation experiments 

Load_GeneletRepressilator_Data.m generates the system identification iddata objects for use in the structured system identification routine.  

DirectQPest_NoPoleZeroCancellations.m executes direct dynamical structure estimation, using an extension of routines from the MATLAB system identification toolbox.  The problem is set up to perform structured identification, via randomized system initializations.  The routine is iterated many times until a best fit is found.  This approach only works on the in vitro transcriptional event detector set, as the perturbation data has a response profile similar to an asymptotically stable linear system.  For unstable data, use Algorithm 1, presented in the paper (also found in the event detector directory). 

gen_pretty_plots.m is a helper function to generate labeled visualizations of time-series traces. 

Gen_Pretty_DSFGraphs.m is an auxiliary function that plots the impulse response functions of each edge in Q(s) as a dynamic graph and saves the animation as a movie in this directory. 
