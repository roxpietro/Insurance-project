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

rates_up = xlsread('EIOPA_RFR_20210331_Term_Structures', 7, 'S11:S30');
%generate the forward
fwd_up = spot2fwd(rates_up,T);

% Bond prices:
Bond_up=BondPricing(fwd_up,T,BT,spread); 
%secondo me va ricalcolato lo spread perchè i rates si alzano anche se dice
%che era costante, bho non ho capito la consegna. 

% The Equity:
Equity_up=EquitySimulation(N,S0,fwd_up,sigma,T);

B0_up=Bond_up(1); 
%non mi torna tanto che il bond adesso valga diverso da 800 cioè capisco
%che deve sempre darmi 1000 però bho è un po' strano
F_up=Bond_up+Equity_up;    % Ft=S_t+B_t

% Case A:
Asset_up=B0_up+S0;      

[Liabilities_A_up,DurL_A_up] = ComputeLiabilities(F_up,rates_up,C0,T,lx,qx,flag_A);

[BOF_A_up,dBOF_A_up,SCR_A_up] = SolvencyComputation(Asset_up,BOF_A_plain,Liabilities_A_up);


% Case B:
 Asset_up=B0_up+S0; 
[Liabilities_B_up,DurL_B_up]=ComputeLiabilities(F_up,rates_up,C0,T,lx,qx,flag_B);    

[BOF_B_up,dBOF_B_up,SCR_B_up]=SolvencyComputation(Asset_up,BOF_B_plain,Liabilities_B_up);

%%  Downward interest rate Scénario :

rates_down = xlsread('EIOPA_RFR_20210331_Term_Structures',8, 'S11:S30');
%generate the forward
fwd_down = spot2fwd(rates_down,T);

% Bond prices
Bond_down=BondPricing(fwd_down,T,BT,spread);

% The Equity
Equity_down = EquitySimulation(N,S0,fwd_down,sigma,T);

B0_down=Bond_down(1);

F_down=Bond_down+Equity_down;   % Ft=S_t+B_t

% Case A:
Asset_down=B0_down+S0;   

[Liabilities_A_down,DurL_A_down]=ComputeLiabilities(F_down,rates_down,C0,T,lx,qx,flag_A);

[BOF_A_down,dBOF_A_down,SCR_A_down]=SolvencyComputation(Asset_down,BOF_A_plain,Liabilities_A_down);

% Case B:
Asset_down=B0_down+S0;

[Liabilities_B_down,DurL_B_down]=ComputeLiabilities(F_down,rates_down,C0,T,lx,qx,flag_B);

[BOF_B_down,dBOF_B_down,SCR_B_down]=SolvencyComputation(Asset_down,BOF_B_plain,Liabilities_B_down); 

%% SCR IR:

% Case A:
SCR_A_interest = max(SCR_A_down,SCR_A_up);
% Case B:
SCR_B_interest = max(SCR_B_down,SCR_B_up);

%% Stressed Equity Scénario :
 
Equity_shock_type1=0.39;

%non mi ricordo cosa avevamo deciso sull'equity adj
S0_shocked = S0*(1-Equity_shock_type1);
 
% Bond prices
Bond=BondPricing(fwd,T,BT,spread);

% The Shocked Equity
Equity_shock =EquitySimulation(N,S0_shocked,fwd,sigma,T);

F_Equity=Bond+Equity_shock;           % Ft=S_t+B_t

% Case A:

Asset_equity=B0+S0_shocked;

[Liabilities_A_Equity,DurL_A_Equity]=ComputeLiabilities(F_Equity,rates,C0,T,lx,qx,flag_A);   

[BOF_A_Equity,dBOF_A_Equity,SCR_A_Equity]=SolvencyComputation(Asset_equity,BOF_A_plain,Liabilities_A_Equity);

% Case B:

Asset_equity=B0+S0_shocked;

[Liabilities_B_Equity,DurL_B_Equity]=ComputeLiabilities(F_Equity,rates,C0,T,lx,qx,flag_B); 

[BOF_B_Equity,dBOF_B_Equity,SCR_B_Equity]=SolvencyComputation(Asset_equity,BOF_B_plain,Liabilities_B_Equity);

%% Stressed spread Scénario :
 
% Credit quality is AAA which is the column 0 and Maturity is 10Y
 
 
f_spread=0.097+0.005*(T-5);                           % stress for duration between 5 and 10 and AAA quality score 
 
B0_shocked=(1-f_spread)*B0;                           % Shocked Value of the zero coupon due to the shock on spread
spread_shocked= -log(B0_shocked/BT) / T - log( 1+rates(end) ); 

Bond_shock=BondPricing(fwd,T,BT,spread_shocked);     % Bond prices:


F_spread=Bond_shock+Equity_plain;                     % Ft=S_t+B_t

% Case A:

Asset_spread=B0_shocked+S0; %cancellerei

[Liabilities_A_spread,DurL_A_spread]=ComputeLiabilities(F_spread,rates,C0,T,lx,qx,flag_A);

Delta_Liabilities_A = Liabilities_A_spread-Liabilities_A_plain;
dBOF_A_spread = f_spread*B0 + Delta_Liabilities_A;
BOF_A_spread = BOF_A_plain-dBOF_A_spread;
SCR_A_spread = max(dBOF_A_spread,0);  

% Case B:
Asset_spread=B0_shocked+S0; %cancellerei

[Liabilities_B_spread,DurL_B_spread]=ComputeLiabilities(F_spread,rates,C0,T,lx,qx,flag_B);

Delta_Liabilities_B=Liabilities_B_spread-Liabilities_B_plain;
dBOF_B_spread=f_spread*B0 +Delta_Liabilities_B;
BOF_B_spread=BOF_B_plain-dBOF_B_spread;
SCR_B_spread=max(dBOF_B_spread,0);

