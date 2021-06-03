% This script executes structured system identification to estimate the dynamical structure function 
% on time-series data generated using Dr. Jongmin Kim's T7 repressilator.  
%  The data was generated using a series of "step input" experiments,
%  conformal with the identifability conditions outlined in Warnick and
%  Goncalves' 2008 paper. 
%
%
% The current approach uses tfest(IDDATA_OBJ) where IDDATA_OBJ is a merged
% iddata object of the different input conditions. 
% Structured parameters for identification are set to either free or
% diagonal based on the 2008 Goncalves and Warnick paper.  An internal
% network structure matrix 
% p x p Q(s) and a control structure matrix P(s) parameterize the network model defined on measured states and inputs.
% Q(s) is assumed to be a square p x p matrix with 0s on the diagonal, to
% ensure the interpretation of causality of the off-diagonal entries.  P(s)
% is a diagonal marix denoting independent input action on each respective output. 
%
%
% Data is processed to ensure that the (Q,P)(s) identified will quantify
% crosstalk in the biological circuit .  

% The identification routine uses an optimum of 5 methods: 1) least-squares
% nonlinear optimization, 2) Levenberg-Marquardt, 3) Gauss-Newton, 4) Subspace Gauss-Newton, and
% 5) gradient descent.  The hyper parameter for this optimization routine
% is the Smith-McMillan Structural Degree of the Dynamical Structure
% Function.   This parameter is varied in a for loop, using a prior
% knowledge of the structure of the distribution of the maximal pole degree
% across Q versus P.  Q is assumed to have uniform maximal pole degree,
% while P is assumed to have distinct uniform maximal pole degree.   These
% assumptions reflect the underlying algebraic structure of Q and P. 
% The algorithm is initialized 100 times with a randomly generated stable state-space model for each pair of
% hyper-parameters on Q and P and for each initialization, optimized to
% provide a best fit structured transfer function model.  The metric for
% optimization used to identify the best hyper-parameter and initialization is the percent fit.  
% 
%


% 
% Script Developed by 
% Enoch H. Yeung, July 2020 
% Center for Control, Dynamical Systems, and Computation
% Department of Mechanical Engineering, University California, Santa
% Barbara
% 

% MATLAB hygiene 
set(0,'DefaultFigureWindowStyle','Normal')
clc 
close all
clear all

%
%
%

DISPLAY_NETWORKSTRUCTURE=true;
DEBUG_SPLASH = true; 
LOAD_DATA_HANDLE = 'Load_GeneletRepressilator_Data';
CIRCUIT_HANDLE = 'Genelet Repressilator';
Ydim = 3; % 3 Outputs in the Repressilator 
Udim = 3; % 3 Inputs in the Repressilator 

if DEBUG_SPLASH
    tic;
end

eval(LOAD_DATA_HANDLE);
 

lbord = 1;
ubord = 1;

opts = tfestOptions();
opts.SearchOption.Tolerance = 1e-16;
opts.SearchMethod = 'auto';
opts.Display = 'on';
opts.Regularization.Lambda=20.0;
opts.SearchOptions.MaxIterations = 100;


OptIndex= Inf*ones(Ydim+Udim,1);
LocalOptVec = OptIndex;
%QPunbiased = cell(ubord+1,ubord+1,ubord+1, ...
%                    ubord+1,ubord+1,ubord+1,...
%                    ubord+1,ubord+1,ubord+1);
%Erro = QPunbiased;
FittoBeat = 0;

for q_poleindex = lbord:ubord
    for p_poleindex = lbord:ubord
        cycle_init=true;
        Max_Tries = 100;
        Num_Tries = 0;
        while (cycle_init == true && Num_Tries < Max_Tries)
                                 Num_Tries = Num_Tries + 1;
                                 if (Ydim~=Udim)
                                     disp("Warning! Your number of inputs should match the number of outputs, see the IEEE Transactions Automatic Control Paper by Dr. Yuan");
                                 end
                                 
                                 %Structured Model Initialization (Random)
                                 InitModel_Raw = zeros(Ydim,Ydim+Udim)*tf('s');  %initialize an empty 0s structure of TFs. 
                                 for row_index = 1:Ydim
                                     for col_index = 1:Ydim+Udim
                                         if col_index<= Ydim && row_index~=col_index  
                                            temp_ss_sys = rss(q_poleindex);
                                            temp_ss_sys.D = temp_ss_sys.D-temp_ss_sys.D;
                                            InitModel_Raw(row_index,col_index) = tf(temp_ss_sys);
                                         end
                                         if col_index> Ydim && (row_index+3)==col_index
                                             temp_ss_sys = rss(p_poleindex);
                                             temp_ss_sys.D = temp_ss_sys.D-temp_ss_sys.D;
                                             InitModel_Raw(row_index,col_index) = tf(temp_ss_sys);
                                         end
                                         InitModel = idtf(InitModel_Raw);
                                         
                                     end
                                 end

                                 for row_index = 1:Ydim
                                        for col_index = 1:Ydim+Udim
                                            if row_index == col_index && col_index <=Ydim
                                                InitModel.Structure(row_index,col_index).Numerator.Value=0.0;
                                                InitModel.Structure(row_index,col_index).Numerator.Free = false; 
                                            else
                                                InitModel.Structure(row_index,col_index).Numerator.Free=1;
                                                InitModel.Structure(row_index,col_index).Denominator.Free(2:end) = 1;

                                            end
                                            if col_index >Ydim 
                                                if (row_index+3) == col_index 
                                                    InitModel.Structure(row_index,col_index).Numerator.Free=1;
                                                    InitModel.Structure(row_index,col_index).Denominator.Free(2:end)=1;
                                                else
                                                    InitModel.Structure(row_index,col_index).Numerator.Value=0.0;
                                                    InitModel.Structure(row_index,col_index).Numerator.Free=false;                
                                                end

                                            end
                                        end
                                 end
                                 try 
                                     tempsys = tfest(totdatobj,InitModel,opts);

                                     %QPunbiased{q12pindex+1,q13pindex+1,p11index+1,q21pindex+1,q23pindex+1,p22index+1,q31pindex+1,q32pindex+1,p33index+1} = tempsys;
                                     ThisErrorObj =resid(tempsys,totdatobj,'corr');
                                     %Erro{q12pindex+1,q13pindex+1,p11index+1,q21pindex+1,q23pindex+1,p22index+1,q31pindex+1,q32pindex+1,p33index+1} = ThisErrorObj;
                                     if min(tempsys.Report.Fit.FitPercent)>FittoBeat  %norm(ThisErrorObj.OutputData,normnum) < ErrotoBeat
                                         %0<min(tempsys.Report.Fit.LossFcn) <ErrotoBeat %
                                         %ErrotoBeat = tempsys.Report.Fit.LossFcn;%norm(ThisErrorObj.OutputData,normnum);
                                         FittoBeat = min(tempsys.Report.Fit.FitPercent);
                                         optsys = tempsys; 
                                         figure(5)
                                         compare(totdatobj,tempsys)
                                         shg
                                         disp(tempsys.Report.Fit.FitPercent)
                                         
                                         
                                         if DISPLAY_NETWORKSTRUCTURE
                                         
                                             Qfin = optsys(1:Ydim,1:Ydim);
                                             Pfin = optsys(1:Ydim,Ydim+1:Ydim+Udim); 

                                             for i = 1:Ydim
                                                 for j = 1:Ydim
                                                     syms s
                                                     [num,den] =  tfdata(Qfin(i,j));
                                                     Q_sym(i,j) =poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
                                                 end
                                             end


                                             digits 32
                                             syms tdum
                                             iQ_nom = ilaplace(Q_sym,tdum);

                                             t = linspace(0,200,300);
                                             for k = 1:numel(t)
                                                 for i = 1:size(iQ_nom,1)
                                                     for j = 1:size(iQ_nom,2)
                                                         iQ_nom_num(i,j,k) = double(subs(iQ_nom(i,j),tdum,t(k)));
                                                     end
                                                 end
                                             end

                                             figure(6)
                                             scale_iQ = iQ_nom_num(:,:,1:end);
                                             

                                             ttemp = t;

                                             for i = 1:size(iQ_nom_num,1)
                                                 for j = 1:size(iQ_nom_num,2)
                                                     h_plot = subplot(3,3,j+3*(i-1));
                                                     h_plot_temp = plot(ttemp(1:end),reshape(scale_iQ(i,j,:),numel(t),1),'r.','MarkerSize',20);
                                                     hold on
                                                     plot(ttemp,zeros(numel(ttemp),1),'-b','LineWidth',3)
                                                     hold off
                                                     
                                                     legend(strcat('Q_{',num2str(i),',',num2str(j),'}'))
                                                     xlim([min(ttemp) max(ttemp)])
                                                     set(gca,'FontSize',20)
                                                     set(gca,'xscale','log')
                                                     ylim([-.1 .1])
                                                 end
                                              end
                                         
                                         end
                                         
                                     end
                                     
                                 catch e %e is an MException struct
                                    fprintf(1,'The identifier was:\n%s',e.identifier);
                                    fprintf(1,'There was an error! The message was:\n%s',e.message);
                                    disp('Failed loss computation')
                                 end

                                 
        end
    end
end

colormat= {[0 0 .8],[0 .8 0],[.8 0 0]};
set(0,'DefaultFigureWindowStyle','docked')

figure(7)
fig_1 = gen_pretty_plots(totdatobj,optsys,tm/60,colormat,2)

figure(5)
compare(totdatobj,optsys,['r.-']);


Qfin = optsys(1:3,1:3);
Pfin = optsys(1:3,4:6); 

for i = 1:3  
    for j = 1:3
        syms s
    [num,den] =  tfdata(Qfin(i,j));
    Q_sym(i,j) =poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
    end
end




temp = vpa(Q_sym,2);
latex(temp)

temp = vpa(ilaplace(Q_sym),2);
latex(temp)
    
digits 32

syms tdum
iQ_nom = ilaplace(Q_sym,tdum);

t = linspace(0,200,300);
for k = 1:numel(t)
    for i = 1:size(iQ_nom,1)
        for j = 1:size(iQ_nom,2)
            iQ_nom_num(i,j,k) = double(subs(iQ_nom(i,j),tdum,t(k)));
        end
    end
end

figure(6)
scale_iQ = iQ_nom_num(:,:,1:end);

ttemp = t;
  
for i = 1:size(iQ_nom_num,1)
    for j = 1:size(iQ_nom_num,2)
        h_plot = subplot(3,3,j+3*(i-1));
        h_plot_temp = plot(ttemp(1:end),reshape(scale_iQ(i,j,:),numel(t),1),'r.','MarkerSize',20)
        hold on 
        plot(ttemp,zeros(numel(ttemp),1),'-b','LineWidth',3)
        hold off
        set(gca,'xscale','log')
        %legend(strcat('Q_{',num2str(i),',',num2str(j),'}'))
        xlim([min(ttemp) max(ttemp)])
        set(gca,'FontSize',20)
        ylim([-.2 .2])

    end
end
save(strcat('7_20_20_direct_structured_est_results',CIRCUIT_HANDLE,'.mat'))
if DEBUG_SPLASH
    toc;
end


colormatrix  = [colormat{1};colormat{2};colormat{3}];
This_Movie_Handle = Gen_Pretty_DSFGraphs(iQ_nom_num(:,:,1:100),colormatrix,1e-2);
save(strcat('7_20_20_direct_structured_est_results',CIRCUIT_HANDLE,'.mat'))
