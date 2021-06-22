function   Bond=BondPricing(r,T,BT,spread)


%function that computes the Bond prices at different instants
%Output: 

% T          : Time to maturity
% N          : Face value
% r          : Interest rates from EIOPA 20200331
% spread     : Spread

t=(1:T)';                                  % time
DF=(1+r).^-t;                              % discount factors
fwd_DF=DF(end)./DF;                        % fwd discounts
B0=BT*DF(end).*exp(-T*spread);   %formula sbagliata           % the zero coupon bond price
Bond = [B0;BT*fwd_DF.*exp(-(T-t)*spread)];  % Bond prices

end 
