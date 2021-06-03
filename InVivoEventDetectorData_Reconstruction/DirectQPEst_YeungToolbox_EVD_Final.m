% This script executes structured system identification (Algorithm 1 in Yeung, Kim et al. 2021) to estimate the dynamical structure function 
% on time-series data generated using Dr. Jongmin Kim's T7 repressilator or the in vivo transcriptional event detector.  
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

%
% Script Developed by 
% Enoch H. Yeung, Dec 2020 
% Center for Control, Dynamical Systems, and Computation
% Department of Mechanical Engineering, University California, Santa
% Barbara
% 
set(0,'DefaultFigureWindowStyle','Normal')
clc 
clear all
%
%
%
clear all
DISPLAY_NETWORKSTRUCTURE=false;
DEBUG_SPLASH = false; 

%LOAD_DATA_HANDLE = 'Load_EventDetector_Data_nM'; % Low
%Inducer Gain Dataset
LOAD_DATA_HANDLE = 'Load_EventDetector_Data_uM';  %High Inducer Gain Dataset
%subtraction not applied in the high induction condition for consistency
%with the low induction condition
CIRCUIT_HANDLE = 'EVENT_DETECTOR';



eval(LOAD_DATA_HANDLE);
close all 


if DEBUG_SPLASH
    tic;
end
Ydim = 2; % 3 Outputs in the Repressilator, automate this 
Udim = 2; % 3 Inputs in the Repressilator 

%ThisOrd =3;
opt_fit_error = Inf;
for ThisOrd = 2:20
    %fit_error = Inf ;

    maxhorizon = 4;
    %this_cond_number = Inf; 

    while 1<maxhorizon  && maxhorizon < length(getexp(totdatobj,'Exp1').Y)-ThisOrd -1
        %display(maxhorizon)
        YLHS_Stack = [];
        YRHS_Stack = []; 
        URHS_Stack = []; 
        Y_Data_Matrix = totdatobj.Y;
        U_Data_Matrix = totdatobj.U;

        %Abstract stacking of data 
        for cell_index = 1:numel(Y_Data_Matrix)
            This_Output_Trace = Y_Data_Matrix{cell_index};
            This_Input_Trace = U_Data_Matrix{cell_index};
            for Column_Index = 1:size(This_Output_Trace,2)  % Iterates over the number of outputs Y1, Y2, ... 
                for Row_Index = 1:(maxhorizon)%step_size:(size(This_Output_Trace,1)-ThisOrd)  % Iterates over the number of timepoints up, less the order of the polynomial 
                    Ntimepoints = size(This_Output_Trace,1); 
                    YLHS_Stack = [YLHS_Stack; This_Output_Trace(Ntimepoints-(Row_Index-1) - ThisOrd:Ntimepoints-(Row_Index-1),Column_Index)'];
                    OtherOutputs_Stack = [];
                    for Other_Outputs = 1:Ydim-1  %There are Ydim-1 other outputs that affect a given Yi in Y=QY + PU  (Q is zero in the diagonal by construction but off-diagonal has free parameters) 
                        if Other_Outputs+Column_Index>Ydim
                            OtherY_Index = mod(Other_Outputs+Column_Index,Ydim);
                        else
                            OtherY_Index = Other_Outputs+Column_Index; 
                        end

                        OtherOutputs_Stack = [OtherOutputs_Stack, This_Output_Trace(Ntimepoints-(Row_Index-1) - (ThisOrd)+1:Ntimepoints-(Row_Index-1),OtherY_Index)'];  %Grab all the other Yjs and concatenate as columns                
                    end
                    YRHS_Stack = [YRHS_Stack; [zeros(1*(Column_Index~=1),(ThisOrd)*(Column_Index-1)*(Ydim-1)), OtherOutputs_Stack,zeros(1*(Ydim~=Column_Index),ThisOrd*(Ydim-Column_Index)*(Ydim-1))]  ];  % Concatenate multiple timepoints as rows e.g. Y[1:4] is one row, the next row would be Y[2:5] for a 2nd order z transform-polynomial 
                    Inputs_Stack = []; 
                    for Input_Index = 1:Udim
                        Inputs_Stack = [Inputs_Stack, (Input_Index==Column_Index)*This_Input_Trace(Ntimepoints-(Row_Index-1)-ThisOrd+1:Ntimepoints-(Row_Index-1),Ydim+Input_Index)']; 
                    end
                    URHS_Stack = [URHS_Stack;Inputs_Stack];
                end
                %Finished iterating over timepoints, now to figure out what to do
                %with different column outputs.  Each column requires a new set of
                %free parameters to parameterize the numerator of entries in Q. 
            end
        end

        X = [-YLHS_Stack(:,2:ThisOrd+1),YRHS_Stack,URHS_Stack];
        B = [YLHS_Stack(:,1)];


        %Calculate the model parameters from Xopt*Theta = Bopt 
        %Theta = inv(Xopt'*Xopt)*Xopt'*Bopt;
        Theta = mldivide(X,B);

        %display(strcat('Structural Fit Error: ',num2str(norm(B- X*Theta)/norm(B))))
        if DEBUG_SPLASH
            toc; 
        end

        % Abstract set up Q and P matrices 
        QNumCell = cell(Ydim,Ydim); 
        QDenCell = cell(Ydim,Ydim); 
        PNumCell = cell(Ydim,Udim);
        PDenCell = cell(Ydim,Udim);
        QPNumCell = cell(Ydim,Ydim+Udim);
        QPDenCell = cell(Ydim,Ydim+Udim);
        DenomCoefficients = [1;Theta(1:ThisOrd)]; 

        for Yindex  = 1:Ydim
            QNumCell{Yindex,Yindex} = zeros(1,ThisOrd);
            QDenCell{Yindex,Yindex} = DenomCoefficients';
            QPNumCell{Yindex,Yindex} = zeros(1,ThisOrd); 
            QPDenCell{Yindex,Yindex} = DenomCoefficients';
        end

        for Rowindex = 1:Ydim 
            for Colindex = 1:Udim 
                QNumCell{Rowindex,Colindex} = zeros(1,ThisOrd);
                PNumCell{Rowindex,Colindex} = zeros(1,ThisOrd);
                PDenCell{Rowindex,Colindex} = DenomCoefficients';  

                QPNumCell{Rowindex,Colindex} = zeros(1,ThisOrd);
                QPNumCell{Rowindex,Colindex+Ydim} = zeros(1,ThisOrd);
                QPDenCell{Rowindex,Colindex} = DenomCoefficients';  
                QPDenCell{Rowindex,Colindex+Ydim} = DenomCoefficients';  
            end
        end

        AllNumeratorCoefficients =cell(Ydim^2-Ydim+Udim,1);
        for Output_Index = 1:(Ydim^2-Ydim+Udim)
            ThisNumeratorCoefficients = (Theta(ThisOrd*(Output_Index)+1:ThisOrd*(Output_Index+1)));
            AllNumeratorCoefficients{Output_Index} = ThisNumeratorCoefficients;
        end

        SkipOffset = 0; 
        for Rowindex = 1:Ydim 
            for Colindex = 1:(Ydim)
                if Rowindex == Colindex 
                    SkipOffset = SkipOffset -1; 
                else
                    QNumCell{Rowindex,Colindex} = AllNumeratorCoefficients{(Rowindex-1)*Ydim +Colindex + SkipOffset}';
                    QPNumCell{Rowindex,Colindex} = AllNumeratorCoefficients{(Rowindex-1)*Ydim +Colindex + SkipOffset}';
                    QDenCell{Rowindex,Colindex} = DenomCoefficients';
                    QPDenCell{Rowindex,Colindex} = DenomCoefficients';
                end
            end
        end
        This_Ts = totdatobj.ts{1}(:);
        Q = d2c(tf(QNumCell,QDenCell,This_Ts),'tustin');
        Qd = (tf(QNumCell,QDenCell,This_Ts));
        for Diagindex = 1:(Udim)
            PNumCell{Diagindex,Diagindex} = AllNumeratorCoefficients{Ydim^2-Ydim +Diagindex}';
            PDenCell{Diagindex,Diagindex} = DenomCoefficients';
            QPNumCell{Diagindex,Diagindex+Ydim} = AllNumeratorCoefficients{Ydim^2-Ydim +Diagindex}';
            QPDenCell{Diagindex,Diagindex+Ydim} = DenomCoefficients';
        end

        QP = d2c(tf(QPNumCell,QPDenCell,This_Ts));

        QPd = (tf(QPNumCell,QPDenCell,This_Ts));
        P = d2c(tf(PNumCell,PDenCell,This_Ts),'tustin');
        Pd = (tf(PNumCell,PDenCell,This_Ts));
        
        
        ErrorVec = zeros(numel(totdatobj.Experiments),1);
        for expid_index = 1:numel(totdatobj.Experiments)
            ThisExpIndex = totdatobj.Experiments{expid_index};
            YPred = lsim(QPd,getexp(totdatobj,ThisExpIndex).u);
            YAct = getexp(totdatobj,ThisExpIndex).y;
            ErrorVec(expid_index) = norm(YPred-YAct,'fro');
        end
            

        if sum(ErrorVec) < opt_fit_error
            opt_fit_error = mean(ErrorVec) ;
            Xopt = X; 
            Bopt = B; 
            opt_ord = ThisOrd;
            opt_horizon = maxhorizon;
            Qfin = Qd; 
            Pfin = Pd; 
            Qcont = Q;
            Pcont = P; 
            QPfin = QPd; 
            QPfincont = QP;
            figure(2)
            hold on
            plot(YPred(:,1),'g-')
            plot(YPred(:,2),'r-')
            plot(getexp(totdatobj,'Exp2').y(:,1),'g.')
            plot(getexp(totdatobj,'Exp2').y(:,2),'r.')
            hold off 
        end
        maxhorizon= maxhorizon+1; 
    end


end

for i = 1:Ydim  
    for j = 1:Ydim
        syms s
        [num,den] =  tfdata(Qfin(i,j));
        Q_sym(i,j) =poly2sym(cell2mat(num),s)/poly2sym(cell2mat(den),s);
    end
end


temp = vpa(Q_sym,2);
latex(temp)

temp = vpa(iztrans(Q_sym),2);

latex(temp)

digits 32

syms tdum

iQ_nom = ilaplace(Q_sym,tdum);

t = linspace(0.0,45,600);
for k = 1:numel(t)
    for i = 1:size(iQ_nom,1)
        for j = 1:size(iQ_nom,2)
            iQ_nom_num(i,j,k) = double(subs(iQ_nom(i,j),tdum,t(k)));
        end
    end
end



figure(1)
scale_iQ = iQ_nom_num(:,:,1:end);

ttemp = t;

for i = 1:size(iQ_nom_num,1)
    for j = 1:size(iQ_nom_num,2)
        h_plot = subplot(Ydim,Ydim,j+Ydim*(i-1));
        h_plot_temp = plot(ttemp(1:end),reshape(scale_iQ(i,j,:),numel(t),1),'r.-','MarkerSize',20)
        hold on 
        plot(ttemp,zeros(numel(ttemp),1),'-b','LineWidth',3)
        hold off
        %set(gca,'xscale','log')
        %legend(strcat('Q_{',num2str(i),',',num2str(j),'}'))
        xlim([0.0 6])
        set(gca,'FontSize',20)
        %ylim([-1.0 1.0])
    end
end
% 

% 
PLOT_DYNAMIC_GRAPH = 1
if PLOT_DYNAMIC_GRAPH 
    colormat= {[.5 .5 0],[.8 0 0],[.1 .1 .1]}
    colormatrix  = [colormat{1};colormat{2}]%;colormat{3}];
    This_Movie_Handle = Gen_Pretty_DSFGraphs(iQ_nom_num(:,:,:),colormatrix,{'',''}',1e-2,CIRCUIT_HANDLE);
end


figure(2)
YPred = lsim(QPd,getexp(totdatobj,'Exp2').u)
hold on
plot(YPred(:,1),'g-')
plot(YPred(:,2),'r-')
plot(getexp(totdatobj,'Exp2').y(:,1),'g.')
plot(getexp(totdatobj,'Exp2').y(:,2),'r.')





for expid_index = 1:numel(totdatobj.Experiments)
    ThisExpIndex = totdatobj.Experiments{expid_index};
    YPred = lsim(QPfin,getexp(totdatobj,ThisExpIndex).u);
    YAct = getexp(totdatobj,ThisExpIndex).y;
    figure();
    hold on;
    plot(YPred(:,1),'g-'); 
    plot(YPred(:,2),'r-')
    plot(getexp(totdatobj,ThisExpIndex).y(:,1),'g.')
    plot(getexp(totdatobj,ThisExpIndex).y(:,2),'r.')
end
