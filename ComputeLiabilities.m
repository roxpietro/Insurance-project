function    [Value,Mac_Duration]=ComputeLiabilities(F,r,C0,T,lx,qx,flag)


% Function that Computes the Liabilities  and the Macaulay Duration 

% time:
t=(1:T)';                
% discount factors:
DF=(1+r).^-t; 

% Compute the fees:
fees = F(1:end-1)*3/100;     
F_prime = F(2:end)-fees;     


if flag == "CaseA"
    C = max(C0,F_prime)
else
    C = F_prime
end

Value=0;
Mac=0;
P=1;    
Pi=[];  
    
for i=1:T 
   %Probability of paying without taking into consideration staying alive at maturity
   Pi = [Pi,(qx(i)+lx(i)*(1-qx(i)))*P]; 
      
   % No lapse No dying in the trajectory
   P=P*(1-qx(i))*(1-lx(i));     
      
   % We pay in case of death/lapse and survivor at maturity, 
   Value = Value + DF(i)*Pi(i)*C(i) + (i==T) * DF(T)*P*C(T); 
      
   % Computing the sum(time*DF*CF*Prob)
   Mac = Mac + t(i)*DF(i)*Pi(i)*C(i) + (i==T) * t(T)*DF(T)*P*C(i); 
end

%  Macaulay Duration
Mac_Duration = Mac/Value;  

end 
