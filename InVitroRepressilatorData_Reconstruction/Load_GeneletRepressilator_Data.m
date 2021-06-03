
load('InVitroRepressilatorData.mat')

normnum = 2; % use L1 norm which weights more equally small deviations.

%The experiment for 08-25-2012 had 3 cuvettes corresponding to 3 different enzyme concentrations (P and H). 
%The cuvette 1 had the most promising response curves.  All cuvettes
%featured a series of 3 perturbations, one for each input. 
%Notation for variables: 
%---------------------------
%c1 is cuvette 1, y_{ij}c_k corresponds to output i, input j (associated with time-window j) and cuvette k. 
%----------------------------



% Begin Data Normalization and Zero-Appending :   
% - -  - - - - - - - - - - -  - - - - - - - 
% - -  - - - - - - - - - - -  - - - - - - - 
% Subtracts out everything up to the input stimulus time,
% thus correcting for any drift from the channels. Such drift is generally assumed to
% be background subtracted in linear time-invariant (LTI) transfer function
% models.  LTI identification routines also require a  sequence of
% non-perturbed, centered trajectories in order to function properly.  The
% premise is that the system is at steady-state, so we append these here.
% A biological network operating in a closed system will experience
% subspace attractors, but the whole system experiences drift that must be
% corrected via the presence of control samples or subtraction of the initial condition.  
% - -  - - - - - - - - - - -  - - - - - - - 
% - -  - - - - - - - - - - -  - - - - - - - 


%  --- Define Window Parameters ---

u1it = 114; %u1 input time 
u2it = 155; %u2 input time 
u3it = 213; %u3 input time 
u4it = 320; %estimated u4 input time (this input was an experimental input that Dr. Kim and Dr. Yeung introduced into the T7 reaction system long after critical data was collected)
%  u4it data is not actually required for estimation of Q and P of the repressilator . 
u1start = u1it; %100 is approximate start time. 
u2start = u2it;  % 160.1 is approximate start time, 165 is start time w/o NMP behavior 
u3start = u3it;  % 230.3 is approximate start time, 240 is start time w/o NMP behavior. 
horlength1 = u2it-u1it;  %length of time-horizon that coincides with u1 acting in isolation 
horlength2 = u3it-u2it;   % length of time-horizon that coincides with u2 acting in isolation 
horlength3 = u4it-u3it;   % length of time-horizon that coincides with u3 acting in isolation 
pren =100; %100 for est. 
nmp = 0;

prez = zeros(1,pren);% Define a length of 
pret = linspace(1,pren,pren);

%Horizon specific smoothing; this ensures we don't accidentally smooth the transients
%associated with steps.  Notall data needed to be smoothed 

sw = 10;  % Smooth/Averaging Window 
y1c1(120:213) = smooth(y1c1(u1it+6:u3it),sw)';
y1c1(230:253) = smooth(y1c1(u3it+17:253),sw)';
y1c1(254:320) = smooth(y1c1(254:320),sw)';


y2c1(240:360) = smooth(y2c1(240:360),sw)';

y3c1(110:155) = smooth(y3c1(110:155),sw)';
y3c1(180:220) = smooth(y3c1(180:220),sw)';
y3c1(229:500) = smooth(y3c1(229:500),sw)';


% - - - - - - - - - - -  - - - - - - - - - 
% Define all input 1 drift corrected outputs (and concomitant time stamp
% vectors)
% - - - - - - - - - - -  - - - - - - - - - 
y11c1 = [prez, y1c1(u1start:u1start+horlength1-1)-[y1c1(u1start:u1start-1), y1c1(u1start)*ones(1,horlength1)]];
y21c1 = [prez, y2c1(u1start:u1start+horlength1-1)-[y2c1(u1start:u1start-1), y2c1(u1start)*ones(1,horlength1)]];
y31c1 = [prez, y3c1(u1start:u1start+horlength1-1)-[y3c1(u1start:u1start-1), y3c1(u1start)*ones(1,horlength1)]];

t11c1 = [pret, pren+t1c1(u1start:u1start+horlength1-1)];
t21c1 = [pret, pren+t2c1(u1start:u1start+horlength1-1)];
t31c1 = [pret, pren+t3c1(u1start:u1start+horlength1-1)];

% - - - - - - - - - - -  - - - - - - - - - 
% Define all input 2 drift corrected outputs (and concomitant time stamp
% vectors)
% - - - - - - - - - - -  - - - - - - - - - 

y12c1 = y1c1(u2start:u2start+horlength2-1)-y1c1(u1start)*ones(1,horlength2);
y22c1 = y2c1(u2start:u2start+horlength2-1)-y2c1(u1start)*ones(1,horlength2); 
y32c1 = y3c1(u2start:u2start+horlength2-1)-y3c1(u1start)*ones(1,horlength2); 

t12c1 = t1c1(u2start:u2start+horlength2-1);
t22c1 = t2c1(u2start:u2start+horlength2-1);
t32c1 = t3c1(u2start:u2start+horlength2-1);

% - - - - - - - - - - -  - - - - - - - - - 
% Define all input 3 drift corrected outputs (and concomitant time stamp
% vectors)
% - - - - - - - - - - -  - - - - - - - - - 

y13c1 = y1c1(u3start:u3start+horlength3-1)-y1c1(u1start)*ones(1,horlength3);
y23c1 = y2c1(u3start:u3start+horlength3-1)-y2c1(u1start)*ones(1,horlength3);
y33c1 = y3c1(u3start:u3start+horlength3-1)-y3c1(u1start)*ones(1,horlength3);

t13c1 = t1c1(u3start:u3start+horlength3-1);
t23c1 = t2c1(u3start:u3start+horlength3-1);
t33c1 = t3c1(u3start:u3start+horlength3-1);



y1c1m = [y11c1, y12c1, y13c1]; % m stands for a merged dataframe (Output 1, Switch T_{23}). 
y2c1m = [y21c1, y22c1, y23c1];% m stands for a merged dataframe (Output 2, Switch T_{12}).
y3c1m = [y31c1, y32c1, y33c1];% m stands for a merged dataframe (Output 3, Switch T_{31}). 
tm = [t11c1 t12c1 t13c1];

u11=10.8*.16*[prez, ones(1,horlength1+horlength2+horlength3 ) ]; % 0.16 µL, added 10.8 uM U3 RNA Inhibitor 3-Analogue
u22=0.16*10.7*[prez,  zeros(1,horlength1)       ,ones(1,horlength2+horlength3)]; %0.16 uL, added 10.7 uM U2 RNA Inhibitor 2-Analogue
u33=0.2*10.5*[prez,  zeros(1,horlength1+horlength2),     ones(1,horlength3)]; %added 0.2 uL, added 10.5 uM U1 RNA Inhibitor 1-Analogue

%Instantitate a IdData Object with All Data in Preparation for System ID
totdatobj = iddata([y1c1m',y2c1m',y3c1m'],[y1c1m',y2c1m',y3c1m',u11' u22' u33']); 

set(0,'DefaultFigureWindowStyle','docked')
figure(1)
plot(totdatobj(:,1,:),'b.--')
figure(2)
plot(totdatobj(:,2,:),'g.--')
figure(3)
plot(totdatobj(:,3,:),'r.--')
figure(4) 
h_plot = plot(totdatobj(:,[1],1),'b.--',totdatobj(:,[2],2),'g.--',totdatobj(:,[3],3),'r.--');


