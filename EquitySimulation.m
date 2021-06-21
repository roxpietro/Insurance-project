function    SimEquity=EquitySimulation(Nbtraj,Nbstep,S0,r,sigma,T)

%function that simulate the equity using the GBM model 

% OUTPUT : Matrix of the simulated path of the equity

% INPUT:
%Nbtraj: Number of trajectories
%Nbstep: Number of step to simulate
%S0 :    Equity at time t=0
%sigma:  Volatility
%r:      rates from the EOPIA 
%T:      Period

%% Generate the B.M:
dt=T/Nbstep;
dW = sqrt(dt)*randn(Nbtraj,Nbstep);  

%% Matrices used to store the approximations
SimEquity=zeros(Nbtraj,Nbstep+1);
SimEquity(:,1) = S0*ones(Nbtraj,1);

%% Generate the fwd rates:
t=(1:T)';                % time
DF=(1+r).^-t;            % discount factors
Zrates=-log(DF)./t;      % zero rates
time_grid=(dt:dt:T)';    % grid of time for simulation
Zrates_grid=interp1(t,Zrates,time_grid,'linear','extrap'); % linear interpolation of the zero rates 
discounts_grid= exp(-time_grid.*Zrates_grid);    
fwd_discount=discounts_grid(2:end)./discounts_grid(1:end-1);  
Rates=[Zrates_grid(1);-log(fwd_discount)/dt];               % DFt+dt = DFt * exp[rt*dt]

%% Compute Equity using GBM

for t=1:Nbstep
SimEquity(:,t+1) = SimEquity(:,t).*exp((Rates(t)-0.5*sigma^2)*dt+sigma*dW(:,t));
end


end