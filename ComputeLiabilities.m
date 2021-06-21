function    [Value,Mac_Duration]=ComputeLiabilities(F,r,C0,T,lx,qx,flag)


% Function that Computes the Liabilities  and the Mecaulay Duration 


% time:
t=(1:T)';                
% discount factors:
DF=(1+r).^-t;            
% Compute the deducted fees:
fees=F(1:end-1,:)*1.5/100;     %we have supposed that the fees are taken on a yearly basis
F_prime = F(2:end,:)-fees;     %F_t'=F_t-f_t where the fees f_t=F_(t-1)*1.5%



if flag== "CaseA"
    Value=0;
    Mac=0;
    C=max(C0,F_prime);
    P=1;    
    Pi=[];  
    for i=1:T
      %Probability of paying without taking into consideration staying alive at maturity
      Pi=[Pi,(qx(i)+lx(i)*(1-qx(i)))*P]; 
      % No lapse No dying in the trajectory
      P=P*(1-qx(i))*(1-lx(i));     
      % we pay in case of death/lapse and survivor at maturity, 
      Value=Value+DF(i)*Pi(i)*mean(C(i,:))+DF(T)*P*mean(C(T,:))*(i==T); 
      % Computing the sum(time*DF*CF*Prob)
      Mac=Mac+t(i)*DF(i)*Pi(i)*mean(C(i,:))+t(T)*DF(T)*P*mean(C(i,:))*(i==T); 
    end
      Mac_Duration=Mac/Value;  %  Mecaulay Duration

end 


if flag== "CaseB"
    Value=0;
    Mac=0;
    C=F_prime;
    P=1;    
    Pi=[];  
    for i=1:T
     % We pay if the insured die or lapses and with no lapsing/no dying for the past instants
     Pi=[Pi,(qx(i)+lx(i)*(1-qx(i)))*P]; 
     % No lapse No dying in the trajectory
     P=P*(1-qx(i))*(1-lx(i));     
     % we pay in case of death/lapse and survivor at maturity 
     Value=Value+DF(i)*Pi(i)*mean(C(i,:))+DF(T)*P*mean(C(T,:))*(i==T); 
     % Computing the sum(time*DF*CF*Prob)
     Mac=Mac+t(i)*DF(i)*Pi(i)*mean(C(i,:))+t(T)*DF(T)*P*mean(C(i,:))*(i==T); 
    end
     Mac_Duration=Mac/Value;    %  Mecaulay Duration
end

end 