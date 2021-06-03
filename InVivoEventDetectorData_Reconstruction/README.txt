### Subrepo manifest: 

*.mat files contain experimental data from the in vitro repressilator validation experiments 

Load_EventDetector_Data_##.m generates the system identification iddata objects for use in the structured system identification routine.  

DirectQPest_YeungToolbox_EVD_Final.m executes direct dynamical structure estimation, using Algorithm 1 from the biorXiv paper: https://www.biorxiv.org/content/10.1101/2021.03.10.434835v2. 

The data for these event detector dynamics resembles an unstable system, since the response of YFP and RFP genes during the relevant time window occurs over exponential growth phase.  The presence or absence of repression in a given network component thus must be deduced from the rates of protein production, rather than by examining the convergence to a steady state, as convergence to a steady-state fluorescence value occurs past log-phase, when cells are dying in culture.   

gen_pretty_plots.m is a helper function to generate labeled visualizations of time-series traces. 

Gen_Pretty_DSFGraphs.m is an auxiliary function that plots the impulse response functions of each edge in Q(s) as a dynamic graph and saves the animation as a movie in this directory. 
