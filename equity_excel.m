function    SimEquity=equity_excel(N,S0,fwd,sigma,T)

%function that simulate the equity using the GBM model 

% OUTPUT : vector of the simulated equity per year

% INPUT:
%N:      Number of samples
%S0:     Equity at time t=0
%sigma:  Volatility
%fwd:    forward from the EOPIA 
%T:      Period

%% Generate the B.M


%% Inizialitation
SimEquity=zeros(T+1,1);
SimEquity(1)=S0;

%% Compute Equity using GBM

for t=2:T+1
SimEquity(t) = SimEquity(t-1).*exp( (fwd(t-1)-0.5*sigma^2));
end

end