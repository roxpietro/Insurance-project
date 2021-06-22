function fwd = spot2fwd(r,T)

t=linspace(1,T,T)'; % time every year  
fwd =-log( ( 1+r(1:end-1) ).^(t(1:end-1) ) ./ ( 1+r(2:end) ).^(t(2:end)) ); %forward rates
fwd=[ -log( (1+r(1))^(-1) ); fwd];

end
