%GROUP 7 SII project

close all;
clear all;
clc;
%%  Data

C0=1000;    % Insured capital 
F0=1000;    % The value of the fund at t_0
B0=800;     % Zero coupon bond price 
BT=1000;    % Face amount 
T=20;       % Maturity
S0= 200;    % Equity price at t_0
sigma=0.15; % Volatility 

%%
% Rates from EIOPA IT with VA 31.03.21
rates = xlsread('EIOPA_RFR_20210331_Term_Structures',4,'S11:S30'); 

%generate the forward
fwd = spot2fwd(rates,T);

% Probability of death (per thousand) ISTA 2017
qx=xlsread('ISTAT 2018 male',1,'E68:E87')/1000;

% Flat annual lapse rates 
lx=0.05*ones(size(qx)); 

%% BASIC CASE
%% Equity

N=100000;    %Number of samples

Equity_plain=EquitySimulation(N,S0,fwd,sigma,T);

%% The Bond :

spread= -log(B0/BT) / T - log( 1+rates(end) ); 
% Bond prices
Bond_plain=BondPricing(fwd,T,BT,spread);


%% Liabilities & Own fund :

flag_A="CaseA";
flag_B="CaseB";

F0=B0+S0;                           % Asset
F_plain=Bond_plain+Equity_plain;    % Ft=S_t+B_t

% Case A:
[Liabilities_A_plain,DurL_A_plain] = ComputeLiabilities(F_plain,rates,C0,T,lx,qx,flag_A);              % Liabilities 
BOF_A_plain = F0-Liabilities_A_plain;                                                                    % Own fund

% Case B:
[Liabilities_B_plain,DurL_B_plain] = ComputeLiabilities(F_plain,rates,C0,T,lx,qx,flag_B);            % Liabilities
BOF_B_plain = F0-Liabilities_B_plain;                                                                  % Own fund


%% Stressed Scénarios :
%% Upward interest rate Scénario:

rates_up = xlsread('EIOPA_RFR_20210331_Term_Structures', 7, 'S11:S20');
% Bond prices:
Bond_up=BondPricing(rates_up,T,N,spread);
% The Equity:
SimEquity=EquitySimulation(Nbtraj,Nbstep,S0,rates_up,sigma,T);
% Equity on time steps
Equity_up=SimEquity(:,1:1/dt:Nbstep+1)';

B0_up=Bond_up(1);            
F_up=Bond_up+Equity_up;    % Ft=S_t+B_t

% Case A:
Asset_up=B0_up+S0;        % Asset
[Liabilities_A_up,DurL_A_up]=ComputeLiabilities(F_up,rates_up,C0,T,lx,qx,flag_A);         % Liabilities/Duration
[BOF_A_up,dBOF_A_up,SCR_A_up]=SolvencyComputation(Asset_up,BOF_A_plain,F_up,rates_up,C0,T,lx,qx,flag_A);

% Case B:
 Asset_up=B0_up+S0;        % Asset
[Liabilities_B_up,DurL_B_up]=ComputeLiabilities(F_up,rates_up,C0,T,lx,qx,flag_B) ;             % Liabilities/Duration
[BOF_B_up,dBOF_B_up,SCR_B_up]=SolvencyComputation(Asset_up,BOF_B_plain,F_up,rates_up,C0,T,lx,qx,flag_B);




