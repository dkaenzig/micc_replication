function tra = SolveHouseholdTransition(p,n,ss,sst,irfcc)

%%% solves for transition of the value function

% set external parameters
ExternalEconomicParameters

% pre-allocate numerical arrays
dVf = zeros(n.NK,1);
dVb = zeros(n.NK,1);
E   = zeros(n.NK,1);

% pre-allocate equilibrium arrays
tra.V     = zeros(n.NK,p.Nt) ;
tra.E     = zeros(n.NK,p.Nt) ;
tra.sav   = zeros(n.NK,p.Nt) ;
tra.Y     = zeros(n.NK,p.Nt) ;
tra.dispY = zeros(n.NK,p.Nt) ;

% assign terminal value
tra.V(:,end)   = sst.V ;
tra.E(:,end)   = sst.E ;
tra.sav(:,end) = sst.sav ;


% iterate backwards
t = p.Nt-1 ;
while t > 0

    % current disposable income
    Y              = irfcc.Z(t) * n.K.^p.alpha ;
    dispY          = max( Y - irfcc.Delta(t) * n.K , 0.0001 ) ;
    tra.Y(:,t)     = Y ;
    tra.dispY(:,t) = Y - irfcc.Delta(t) * n.K ;
    
    % tomorrow's value
    V = tra.V(:,t+1) ;
    
    % forward difference
    dVf(1:n.NK-1) = ( V(2:n.NK) - V(1:n.NK-1) ) ./ n.dK(1:n.NK-1) ;
    dVf(n.NK,:)   = dispY(n.NK).^(-p.gamma); % will never be used, but impose state constraint a<=amax just in case
    
    % backward difference
    dVb(2:n.NK) = ( V(2:n.NK) - V(1:n.NK-1) ) ./ n.dK(1:n.NK-1) ;
    dVb(1,:)    = dispY(1).^(-p.gamma); % state constraint boundary condition, for stabilitiy
    
    % consumption and savings with forward difference
    Ef = up1(dVf) ;
    sf = dispY - Ef ;
    
    % consumption and savings with backward difference
    Eb = up1(dVb);
    sb = dispY - Eb ;
    
    % consumption and derivative of value function at steady state
    E0 = dispY ;
    
    % dV_upwind makes a choice of forward or backward differences based on
    % the sign of the drift    
    If = sf > 0; % positive drift --> forward difference
    Ib = sb < 0; % negative drift --> backward difference
    I0 = (1 - If - Ib ); % at steady state
    
    E = Ef .* If + Eb .* Ib + E0 .* I0 ;
    tra.E(:,t)   = E ;
    tra.sav(:,t) = dispY - E ;
   
    uflow = u(E) ;

    % construct finite difference matrix for assets
    Sb = - min( sb , 0 ) ./ n.dK; % X
    Sm = - max( sf , 0 ) ./ n.dK + min( sb , 0 ) ./ n.dK ; % Y
    Sf =   max( sf , 0 ) ./ n.dK ; % Z
    
    S = spdiags( Sm(:)                  , 0  , n.NK , n.NK ) ...
      + spdiags( [ Sb(2:n.NK) ; 0 ]    , -1 , n.NK , n.NK ) ...
      + spdiags([ 0 ; Sf(1:n.NK-1) ]   , 1  , n.NK ,n.NK ) ;
    
    L = S  ;
    
    % define matrix to invert
    M = ( 1 / p.dt + p.rho ) * speye(n.NK) - L ;

    % invert linear system
    Vnew =  M \ ( uflow + V / p.dt  ) ;
    
    % assigne new value
    tra.V(:,t) = Vnew ;
    t = t-1 ;

end

% assign outcomes
sst.V   = V ;
sst.E   = E ;
sst.sf  = sf ;
sst.sb  = sb ;
sst.sav = dispY - E ;


end

