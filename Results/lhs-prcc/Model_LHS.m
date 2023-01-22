%% The results should be compared to the PRCC results section in
%% Supplementary Material D and Table D.1 for different N (specified by
%% "runs" in the script below
%
% Follow this: http://malthus.micro.med.umich.edu/lab/usadata/
%
%
clear all;
close all;

%% Sample size N
runs=100;

%% LHS MATRIX  %%
% variable time_points is set in this function call (there is a m file in the directory with this name).
Parameter_settings_LHS;

s_LHS			=LHS_Call(1e-2, s, 50, 0 ,runs,'unif')			; % baseline = 10
% size(s_LHS)
muT_LHS		=LHS_Call(1e-4, muT, 0.2, 0 ,runs,'unif')		; % baseline = 2e-2
% size(muT_LHS)
r_LHS			=LHS_Call(1e-3, r, 50, 0, runs,'unif')			; % baseline = 3e-2
% size(r_LHS)
k1_LHS			=LHS_Call(1e-7,k1,1e-3, 0 ,runs,'unif')			; % baseline = 2.4e-5
% size(k1_LHS)
k2_LHS			=LHS_Call(1e-5, k2, 1e-2, 0, runs,'unif')			; % baseline = 3e-3
% size(k2_LHS)
mub_LHS		= LHS_Call(1e-1 , mub , 0.4 , 0 , runs , 'unif')	; % baseline = 0.24
% size(mub_LHS)
N_LHS			=LHS_Call(1,N,2e3, 0 ,runs,'unif')				; % baseline = 1200
% size(N_LHS)
muV_LHS		=LHS_Call(1e-1,muV,1e1, 0 ,runs,'unif')			; % baseline = 2.4
% size(muV_LHS)
dummy_LHS	=LHS_Call(1,1,1e1, 0 ,runs,'unif')				; % dummy parameter
% size(dummy_LHS)

%% LHS MATRIX and PARAMETER LABELS
LHSmatrix=[s_LHS muT_LHS r_LHS k1_LHS k2_LHS mub_LHS N_LHS muV_LHS  dummy_LHS];

% size(LHSmatrix)

for x=1:runs %Run solution x times choosing different values
    f=@ODE_LHS;
%    x
%    LHSmatrix(x,:)
    [t,y]=ode15s(@(t,y)f(t,y,LHSmatrix,x,runs),tspan,y0,[]); 
    A=[t y]; % [time y]
    %% Save the outputs at ALL time points [tspan]
    %T_lhs(:,x)=Anew(:,1);
    %CD4_lhs(:,x)=Anew(:,2);
    %T1_lhs(:,x)=Anew(:,3);
    %T2_lhs(:,x)=Anew(:,4);
    %V_lhs(:,x)=Anew(:,5);
    
    %% Save only the outputs at the time points of interest [time_points]:
    %% MORE EFFICIENT
    T_lhs(:,x)		=	A(time_points+1,1)		;
    CD4_lhs(:,x)	=	A(time_points+1,2)	;
    T1_lhs(:,x)		=	A(time_points+1,3)	;
    T2_lhs(:,x)		=	A(time_points+1,4)	;
    V_lhs(:,x)		=	A(time_points+1,5)	;
end
%% Save the workspace
save Model_LHS.mat;
% CALCULATE PRCC
% I dont know the value of alpha yet. correct alpha is needed to get correct prcc.
% [prcc sign sign_label]=PRCC(LHSmatrix,V_lhs,1:length(time_points),PRCC_var,alpha);
%
%
%size(LHSmatrix)
%size(V_lhs)
%time_points
%PRCC_var
[prcc sign sign_label]=PRCC(LHSmatrix,V_lhs,1:length(time_points),PRCC_var, 0.5);






