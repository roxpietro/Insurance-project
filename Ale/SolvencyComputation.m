function [BOF,dBOF,SCR]=SolvencyComputation(Asset,BOF_plain,F,rates,C0,T,lx,qx,flag)

%da migliorare perchè calcolo le liabilities anche quando sono fuori -->
%forse meglio passare come input 
[Liabilities,~]=ComputeLiabilities(F,rates,C0,T,lx,qx,flag);                  %  Liabilities 
BOF=Asset-Liabilities;                                                     % Own fund
dBOF=BOF_plain-BOF;
SCR=max(dBOF,0);                                                           % Solvency 

end 
