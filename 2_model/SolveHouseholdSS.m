function sst = SolveHouseholdSS(p,n,sst,init)

%%% solves for steady-state value function

% set external parameters
ExternalEconomicParameters

% pre-allocate equilibrium arrays
dVf = zeros(n.NK,1);
dVb = zeros(n.NK,1);
E   = zeros(n.NK,1);

% initial guesses
Y     = sst.Z * n.K.^p.alpha ;
dispY = max( Y - sst.Delta * n.K , 0.0001) ;
if isempty(init)
    v     = u(Y) / p.rho ;
else
    v = init.V ;
end

% initialize loop
error = 1 ;
j = 1 ;

% loop: converge value function
while error > n.crit && j < n.maxit
    
    % store value
    V = v;
    
    % forward difference
    dVf(1:n.NK-1) = ( V(2:n.NK) - V(1:n.NK-1) ) ./ n.dK(1:n.NK-1) ;
    dVf(n.NK,:)   = dispY(n.NK).^(-p.gamma); % will never be used, but impose state constraint a<=amax just in case
    
    % backward difference
    dVb(2:n.NK) = ( V(2:n.NK) - V(1:n.NK-1) ) ./ n.dK(1:n.NK-1) ;
    dVb(1,:)    = dispY(1).^( -p.gamma ); % state constraint boundary condition, not binding but necessary for stability
    
    % consumption and savings with forward difference
    Ef = up1(dVf);
    sf = dispY - Ef ;
    
    % consumption and savings with backward difference
    Eb = up1(dVb);
    sb = dispY - Eb ;
    
    % consumption and derivative of value function at steady state
    E0 = dispY ;
    
    % dV_upwind makes a choice of forward or backward differences based on the sign of the drift    
    If = sf > 0; % positive drift --> forward difference
    Ib = sb < 0; % negative drift --> backward difference
    I0 = (1 - If - Ib ); % at steady state
    
    E = Ef .* If + Eb .* Ib + E0 .* I0 ;
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
    M = ( 1 / n.Delta + p.rho ) * speye(n.NK) - L ;

    % invert linear system
    V =  M \ ( uflow + V / n.Delta  ) ;
    
    Vchange = V - v;
    v = V ;

    dist(j) = max(max(abs(Vchange(:)))) ;
    error = dist(j) ;
    j = j+1 ;

end
% assign outcomes
sst.V   = V ;
sst.E   = E ;
sst.sav = dispY - E ;


end

