function [sst_GWp, irfg_GWp, ...
          sst_SCCp, irfg_SCCp, ss ] = SolveCounterfactuals(optns,ss,p,zetas,init)


%%% solves for model counterfactuals: transitional dynamics and SCC

% re-set baseline parameters for safety
ExternalEconomicParameters ;

% initial steady-state
ss.R = p.rho + p.Delta ;
ss.K = ( p.alpha * p.Z / ss.R )^(1/(1-p.alpha)) ;
ss.Y = p.Z * ss.K^p.alpha ;
ss.E = ss.Y - p.Delta * ss.K ;
ss.I = p.Delta * ss.K ;
ss.Z = p.Z ;
ss.Delta = p.Delta ;

% set grid for steady-state
sst = ss ;
NumericalParameters ;

% assign pre-specified guess if applicable
if nargin==4
    init_ss      = [] ;
    init_sstWp   = [] ;
    init_sstWb   = [] ;
    init_sstSCCp = [] ;
    init_sstSCCb = [] ;
elseif nargin>=5 && isempty(init)
    init_ss      = [] ;
    init_sstWp   = [] ;
    init_sstWb   = [] ;
    init_sstSCCp = [] ;
    init_sstSCCb = [] ;
elseif nargin>=5 && ~ismember('ss',fieldnames(init))
    init_ss      = [] ;
    init_sstWp   = [] ;
    init_sstWb   = [] ;
    init_sstSCCp = [] ;
    init_sstSCCb = [] ;
else
    init_ss      = init.ss ;
    init_sstWp   = init.sstWp ;
end

% compute steady-state value function
ss2 = SolveHouseholdSS(p,n,ss,init_ss) ;
ss.V = ss2.V ;
ss.ss = ss ;


%%% climate warming counterfactual

% define warming scenario
if strcmp(optns.retro,"retro")
    TD = p.TD ;
else
    WarmingScenario
end

% compute counterfactual and assign results
shocks.zetas           = zetas ;
[ sst_GWp , irfg_GWp ] = Counterfactual(p,ss,shocks,TD,init_sstWp) ;
ss.sstWp = sst_GWp ;


%%% SCC counterfactual

% construct temperature path
TCO2          = p.adj * p.CO2size.current ...
              * ( p.xD(1) * ( exp( -p.xD(2)*p.times ) - exp( -p.xD(3)*p.times ) ) + p.xD(4)*( 1 - exp( -p.xD(5) * p.times ) ) ) ;

% assign temperature path
p.TCO2 = TCO2 ;

% compute counterfactual and assign results
shocks.zetas             = zetas ;
[ sst_SCCp , irfg_SCCp ] = Counterfactual(p,ss,shocks,p.TCO2,ss) ;

% assign correct SCC to GW irf to avoid confusion
irfg_GWp.SCC = irfg_SCCp.SCC ;



end