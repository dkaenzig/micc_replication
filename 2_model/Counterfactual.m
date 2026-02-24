function [sst,irfg] = Counterfactual(p,ss,shocks,TD,init)

% set external parameters for safety
ExternalEconomicParameters

%%% construct cumulative shocks

% initialize
zhatc     = zeros( p.Nt , 1) ;
Deltahatc = zeros( p.Nt , 1) ;

for t=1:p.Nt
    
    zhatct = 0 ;
    for s=1:(t-1)
        zhatct     = zhatct     + p.dt * TD(s) * shocks.zetas(t-s) ;
    end
    zhatc(t)     = zhatct ;

end
irfcc.zhat     = zhatc ;

% transform log changes to levels according to model
irfcc.Z     = p.Z * exp( irfcc.zhat ) ;
irfcc.Delta = ss.Delta * ones(p.Nt,1) ;


%%% terminal steady-state

% terminal steady-state shocks
sst.Z     = irfcc.Z(end) ;
sst.Delta = irfcc.Delta(end) ;

% equilibrium objects
sst.R = p.rho + sst.Delta ;
sst.K = ( p.alpha * sst.Z / sst.R )^(1/(1-p.alpha)) ;
sst.Y = sst.Z * sst.K^p.alpha ;
sst.E = sst.Y - sst.Delta * sst.K ;
sst.Ess = sst.E ; % to avoid overwriting below
sst.I = sst.Delta * sst.K ;

% define grids
NumericalParameters

% value function
if nargin == 4 
    init = [] ;
end
sst = SolveHouseholdSS(p,n,sst,init) ;

% find capital index of initial steady-state capital
Viss      = u(ss.E)/p.rho ;
ss.Viss   = Viss ;

%%% solve for the transition

% solve transition
tra      = SolveHouseholdTransition(p,n,ss,sst,irfcc) ;
irfcc.K0 = ss.K ;
irfg     = IRF(p,n,ss,irfcc,tra) ;

% assign outcomes
irfg.TD         = TD ;
irfg.times      = p.times ;
irfg.zhat       = zhatc ;
irfg.Deltahat   = Deltahatc ;

% "SCC": stock valuation of counterfactual
% only valid as SCC when using the temperature-CO2 pulse
k = FindZeroInterp( ss.V - irfg.V(1) , n.K ) ;
irfg.SCC = ( ss.K - k ) / ss.E * p.WorldCons / ( p.CO2size.current * p.Gt100 ) ;

end