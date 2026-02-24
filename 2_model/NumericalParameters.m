

%%% numerical parameters

% capital grid
n.NK   = 1000 ; 

n.Kmin = 0.0001 ;
n.Kmax = n.Kmin + 1.2 * max( sst.K , ss.K ) ;

n.K  = linspace(n.Kmin , n.Kmax , n.NK )';
n.dK = [ n.K(2:end) - n.K(1:end-1) ; n.K(end) - n.K(end-1) ];

% convergence and smoothing criteria
n.Delta = 1 ; 
n.maxit = 2000;
n.crit  = 1e-8 ; 

% maximum number of steady-state iterations
n.Ir   = 1000 ;


