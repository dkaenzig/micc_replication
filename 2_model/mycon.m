function [c,ceq] = mycon(x,q)

%%% nonlinear constraint that returns the ratio of damages at date T to the maximum of the damage function

T = 10 ;

times = 0:1:299 ;
zeta = exp( -x(2) * times ) .* ( 1 - exp( -x(3) * times ) ) ;
mzeta = max(zeta) ;
 
c = exp( -x(2) * T ) .* ( 1 - exp( -x(3) * T ) ) / mzeta - q ;
ceq = 0 ;


end