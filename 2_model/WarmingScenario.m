
%%% defines warming scenario

Tasymptote  = 1.1 * p.T2100 ;
InitYear    = 2024 ;
phiTD       = - 1 / (2100-InitYear) * log(1-p.T2100/Tasymptote) ;
TD          = Tasymptote * ( 1 - exp( - phiTD * p.times ) ) ;