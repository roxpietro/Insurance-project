function [BOF,dBOF,SCR]=SolvencyComputation(Asset,BOF_plain,Liabilities)


BOF=Asset-Liabilities;                                                     
dBOF=BOF_plain-BOF;
SCR=max(dBOF,0);                                                           

end 
