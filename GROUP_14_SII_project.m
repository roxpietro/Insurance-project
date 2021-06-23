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

% Bond prices
Bond_up=BondPricing(fwd_up,T,BT,spread);  

% The Equity:
Equity_up=EquitySimulation(N,S0,fwd_up,sigma,T);

B0_up=Bond_up(1); 
F_up=Bond_up+Equity_up;

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

F_down=Bond_down+Equity_down;

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

S0_shocked = S0*(1-Equity_shock_type1);

% The Shocked Equity
Equity_shock =EquitySimulation(N,S0_shocked,fwd,sigma,T);

F_Equity=Bond_plain+Equity_shock;

% Case A:
F0_shock=B0+S0_shocked;

[Liabilities_A_Equity,DurL_A_Equity]=ComputeLiabilities(F_Equity,rates,C0,T,lx,qx,flag_A);   

[BOF_A_Equity,dBOF_A_Equity,SCR_A_Equity]=SolvencyComputation(F0_shock,BOF_A_plain,Liabilities_A_Equity);

% Case B:

F0_shock=B0+S0_shocked;

[Liabilities_B_Equity,DurL_B_Equity]=ComputeLiabilities(F_Equity,rates,C0,T,lx,qx,flag_B); 

[BOF_B_Equity,dBOF_B_Equity,SCR_B_Equity]=SolvencyComputation(F0_shock,BOF_B_plain,Liabilities_B_Equity);

%% Stressed spread Scénario :
 
% Credit quality is AAA which is the column 0 and Maturity is 10Y
 
 
f_spread=0.097+0.005*(T-15);                          % stress for duration between 15 and 20 and AAA quality score 
 
B0_shocked=(1-f_spread)*B0;                           % Shocked Value of the zero coupon due to the shock on spread
spread_shocked= -log(B0_shocked/BT) / T - log( 1+rates(end) ); 

Bond_shock=BondPricing(fwd,T,BT,spread_shocked);     % Bond prices:


F_spread=Bond_shock+Equity_plain;
Asset_spread=B0_shocked+S0;

% Case A:

[Liabilities_A_spread,DurL_A_spread]=ComputeLiabilities(F_spread,rates,C0,T,lx,qx,flag_A);

[BOF_A_spread,dBOF_A_spread,SCR_A_spread]=SolvencyComputation(Asset_spread,BOF_A_plain,Liabilities_A_spread); 

% Case B:

[Liabilities_B_spread,DurL_B_spread]=ComputeLiabilities(F_spread,rates,C0,T,lx,qx,flag_B);

[BOF_B_spread,dBOF_B_spread,SCR_B_spread]=SolvencyComputation(Asset_spread,BOF_B_plain,Liabilities_B_spread); 

%% Stressed mortality Scénario

qx_shocked=1.15*qx;

% Case A:
[Liabilities_A_mortality,DurL_A_mortality]=ComputeLiabilities(F_plain,rates,C0,T,lx,qx_shocked,flag_A);

[BOF_A_mortality,dBOF_A_mortality,SCR_A_mortality]=SolvencyComputation(F0,BOF_A_plain,Liabilities_A_mortality);

% Case B:
[Liabilities_B_mortality,DurL_B_mortality]=ComputeLiabilities(F_plain,rates,C0,T,lx,qx_shocked,flag_B);

[BOF_B_mortality,dBOF_B_mortality,SCR_B_mortality]=SolvencyComputation(F0,BOF_B_plain,Liabilities_B_mortality);


%% Stressed lapse Scénario:

lx_up=min(1.5*lx,1);
lx_down=max(0.5*lx,lx-0.2);
lx_mass=[0.4;lx(2:end)];

% Case A:

% lapse up
[Liabilities_A_lapse_up,DurL_A_lapse_up]=ComputeLiabilities(F_plain,rates,C0,T,lx_up,qx,flag_A);

[BOF_A_lapse_up,dBOF_A_lapse_up,SCR_A_lapse_up]=SolvencyComputation(F0,BOF_A_plain,Liabilities_A_lapse_up);

% lapse down
[Liabilities_A_lapse_down,DurL_A_lapse_down]=ComputeLiabilities(F_plain,rates,C0,T,lx_down,qx,flag_A);

[BOF_A_lapse_down,dBOF_A_lapse_down,SCR_A_lapse_down]=SolvencyComputation(F0,BOF_A_plain,Liabilities_A_lapse_down);

% lapse mass
[Liabilities_A_lapse_mass,DurL_A_lapse_mass]=ComputeLiabilities(F_plain,rates,C0,T,lx_mass,qx,flag_A);

[BOF_A_lapse_mass,dBOF_A_lapse_mass,SCR_A_lapse_mass]=SolvencyComputation(F0,BOF_A_plain,Liabilities_A_lapse_mass);



SCR_A_lapse=max([SCR_A_lapse_up,SCR_A_lapse_down,SCR_A_lapse_mass]);


% Case B:

% lapse up
[Liabilities_B_lapse_up,DurL_B_lapse_up]=ComputeLiabilities(F_plain,rates,C0,T,lx_up,qx,flag_B);

[BOF_B_lapse_up,dBOF_B_lapse_up,SCR_B_lapse_up]=SolvencyComputation(F0,BOF_B_plain,Liabilities_B_lapse_up);

% lapse down
[Liabilities_B_lapse_down,DurL_B_lapse_down]=ComputeLiabilities(F_plain,rates,C0,T,lx_down,qx,flag_B);

[BOF_B_lapse_down,dBOF_B_lapse_down,SCR_B_lapse_down]=SolvencyComputation(F0,BOF_B_plain,Liabilities_B_lapse_down);


% lapse mass
[Liabilities_B_lapse_mass,DurL_B_lapse_mass]=ComputeLiabilities(F_plain,rates,C0,T,lx_mass,qx,flag_B);

[BOF_B_lapse_mass,dBOF_B_lapse_mass,SCR_B_lapse_mass]=SolvencyComputation(F0,BOF_B_plain,Liabilities_B_lapse_mass);



SCR_B_lapse=max([SCR_B_lapse_up,SCR_B_lapse_down,SCR_B_lapse_mass]);


%% CAT Scénario:

qx_cat=[qx(1)+0.0015; qx(2:end)];


% Case A:
[Liabilities_A_cat,DurL_A_cat]=ComputeLiabilities(F_plain,rates,C0,T,lx,qx_cat,flag_A);

[BOF_A_cat,dBOF_A_cat,SCR_A_cat]=SolvencyComputation(F0,BOF_A_plain,Liabilities_A_cat);

% Case B:
[Liabilities_B_cat,DurL_B_cat]=ComputeLiabilities(F_plain,rates,C0,T,lx,qx_cat,flag_B);

[BOF_B_cat,dBOF_B_cat,SCR_B_cat]=SolvencyComputation(F0,BOF_B_plain,Liabilities_B_cat);



%% LIFE Solvency Capital Requirement:

Correlations_life=[1 0 0.25; 0 1 0.25; 0.25 0.25 1];

% Case A:
SCR_A_LIFE=[SCR_A_mortality; SCR_A_lapse; SCR_A_cat];

SCR_A_life=sqrt(SCR_A_LIFE'*Correlations_life*SCR_A_LIFE);

% Case B:

SCR_B_LIFE=[SCR_B_mortality; SCR_B_lapse; SCR_B_cat];

SCR_B_life=sqrt(SCR_B_LIFE'*Correlations_life*SCR_B_LIFE);


%% MARKET Solvency Capital Requirement:


Correlations_interest_up=[1 0 0; 0 1 0.75; 0 0.75 1];
Correlations_interest_down=[1 0.5 0.5; 0.5 1 0.75; 0.5 0.75 1];


% Case A:
SCR_A_MARKET=[SCR_A_interest; SCR_A_Equity; SCR_A_spread];
SCR_A_MARKET_interest_up=sqrt(SCR_A_MARKET'*Correlations_interest_up*SCR_A_MARKET);
SCR_A_MARKET_interest_down=sqrt(SCR_A_MARKET'*Correlations_interest_down*SCR_A_MARKET);

SCR_A_market=max(SCR_A_MARKET_interest_down,SCR_A_MARKET_interest_up);

% Case B:
SCR_B_MARKET=[SCR_B_interest; SCR_B_Equity; SCR_B_spread];
SCR_B_MARKET_interest_up=sqrt(SCR_B_MARKET'*Correlations_interest_up*SCR_B_MARKET);
SCR_B_MARKET_interest_down=sqrt(SCR_B_MARKET'*Correlations_interest_down*SCR_B_MARKET);

SCR_B_market=max(SCR_B_MARKET_interest_down,SCR_B_MARKET_interest_up);


%% BASIC Solvency Capital Requirement

Correlations_BSCR=[1 0.25; 0.25 1];

% Case A:
SCR__A=[SCR_A_market SCR_A_life]';
BSCR_A=sqrt(SCR__A'*Correlations_BSCR*SCR__A);

% Case B:
SCR__B=[SCR_B_market SCR_B_life]';
BSCR_B=sqrt(SCR__B'*Correlations_BSCR*SCR__B);


