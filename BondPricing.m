function   Bond=BondPricing(fwd,T,BT,spread)


%function that computes the Bond prices at different instants
%Output: 
%Bond = vector of bond prices per year

% T          : Time to maturity
% BT         : Face value
% fwd        : forward from EIOPA 20200331
% spread     : Spread

Bond = zeros(T+1,1);
Bond(T+1)=BT;
for i=T:-1:1
   Bond(i) = Bond(i+1) * exp( -(fwd(i) + spread) );
end

end 
