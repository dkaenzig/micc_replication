function [p, ss, fit, ...
          sst_GWp, irfg_GWp, ...
          sst_SCCp, irfg_SCCp ] = EstimateAndSolve(optns,p,data,init)


%%% model steady-state

% re-set external parameters for safety
ExternalEconomicParameters ;

% steady-state aggregates
ss.R = p.rho + p.Delta ;
ss.K = ( p.alpha * p.Z / ss.R )^(1/(1-p.alpha)) ;
ss.Y = p.Z * ss.K^p.alpha ;
ss.E = ss.Y - p.Delta * ss.K ;
ss.I = p.Delta * ss.K ;
ss.Z = p.Z ;
ss.Delta = p.Delta ;



%%% log-linearized transitional dynamics

% transition matrix
ss.L = [ [ p.rho                                  , - 1 ] ; ...
         [ - (1-p.alpha)*ss.E*ss.R/(p.gamma*ss.K) , 0   ] ] ;

% find eigenspaces
[M0,D] = eig(ss.L) ; % L = M0*D*inv(M0) ;
M        = inv(M0) ; % L = inv(M)*D*M ;

% slope of consumption as a function of capital
ss.eK = - M(1,1)/M(1,2) ;

% shock vector for z
ss.Ss     = [ ss.Y ; ss.E * ss.R / p.gamma ] ;

% slope of consumption as a function of shocks
ss.eS = - p.dt * exp(- D(1,1) * p.times ) / M(1,2) * M(1,:) ;

% decay of capital
p.kappa = - ( ss.L(1,1)-ss.eK ) ;
    
% construct J
D1       = D(1,1) ;
D1_kappa = D1 + p.kappa;
J        = exp(-D1 * p.times') ...
         .* (exp(D1_kappa * min(p.times,p.times') ) - 1) / D1_kappa ;

% construct cK
eS0 = - 1 / M(1,2) * M(1,:) ;
cK  = exp(-p.kappa.*p.times ) / ss.K ...
   .* (   p.dt * ss.Ss(1) * exp( p.kappa * p.times' ) .* ( p.times' <= p.times ) ...
        - p.dt * eS0 * ss.Ss * J ) ;
cK(1, :) = 0 ;


    
%%% process temperature, output and capital data

% extract data IRFs up to year 10
Tdata  = data.temperature ;
Ypdata = data.gdp/100 ;
Kpdata = data.capital/100 ;

% fill in temperature post horizon 10
% estimate AR(1) to fill in gap
% only matters for expectations, not realizations

% select temperature observations
Ypdata2 = Ypdata(2:end) ;
Tdata2  = Tdata(2:end) ;

% keep only positive temperature values
Ypdata2 = Ypdata2(Tdata2>0) ;
Tdata2  = Tdata2(Tdata2>0) ; 

% compute AR coeffcient
b = mean( log(Tdata2) ) / mean(1:(size(Ypdata2,1))) ;
T = [  Tdata ; exp( b * p.times((1+size(Ypdata,1)):p.Nt) ) ] ;

% matrix to convert from transitory to persistent; corresponds to convolution operation
L = ConvolutionMatrix(p,T) ;

% if use transitory: 
% convert to effect of transitory shock using Sims (1986)
if optns.persistence == "transitory"

    Ypdata = L(1:11,1:11) \ Ypdata ; 
    T      = L \ T ;
    L      = eye(size(L)) ;

end

% not used in estimation, but used for plots
Yp     = [ Ypdata ; zeros(p.Nt-size(Ypdata,1),1) ] ;
Kp     = [ Kpdata ; zeros(p.Nt-size(Kpdata,1),1) ] ;

% assign to outcomes
fit.Yp = Yp ;
fit.Kp = Kp ;



%%% estimation

% estimation horizon and initial guess
MaxYear = optns.MaxYear ;
x0prod = [ - 6 , 0.5 , 0.01 ]' ;
    
% productivity in capital targeting case
LHSY    = Yp ;

% change initial guess if supplied
if nargin >= 4 && ~isempty(init) 
    if ismember('x0',fields(init))
         x0prod = init.xs ; 
    end
end

% damage function
zeta   = @(x) x(1) * exp( -x(2) * p.times ) .* ( 1 - exp( -x(3) * p.times ) ) ;

% implied productivity and output paths
TFP              = @(x) L * zeta(x) ;
outputPersistent = @(x) TFP(x) + p.alpha * cK * TFP(x) ;

% default weights
weights = ones(p.Nt,1) .* ( p.times <= MaxYear) ;

% objective function     
objzeta = @(x) sum( weights .* ( LHSY - outputPersistent(x) ).^2 ) ;

% optimization options
options = optimoptions('fmincon') ;
options.Display = 'none' ;

% set nonlinear contraint, 
% never binding in baseline, only useful for regional damages
Mycon = @(x) mycon(x, p.q) ;

% initial guessL simple OLS LHS: = x1 * RHS
RHS     = outputPersistent( [1 ; x0prod(2:end) ] ) ; % use shape of point estimate
lhs     = LHSY(p.times <= MaxYear) ;
rhs     = RHS(p.times <= MaxYear) ;
ols     = dataset(lhs,rhs) ;
fitols  = fitlm(ols,'lhs ~ rhs-1') ; % -1 to force constant to be zero
xs(1)   = fitols.Coefficients.Estimate(1) ;
xs      = [xs(1) ; x0prod(2:end) ] ;

% run optimization 
xs = fmincon(objzeta , ... objective function
             xs , ...      starting point
             [],[],[],[],...
             [ p.Amin ; p.Bmin ; p.Cmin ] , [ p.Amax ; p.Bmax ; p.Cmax ] ,...  lower and upper bound on parameters
             Mycon , options ) ; % nonlinear contraints and options

% recover damage function and store
zetas     = zeta(xs) ;

% recover implied productivity shocks from possibly persistent temperature path
zhat     = TFP(xs) ;

% recover implied output and capital path to persistent shock
yhat     = outputPersistent(xs) ;
khat     = cK * zhat ;

% store estimates
fit.xs      = xs ;
fit.zetas   = zetas ;
fit.zhat    = zhat ;
fit.Yp_fit  = 100 * yhat ;
fit.Kp_fit  = 100 * khat ;

% solve for all household counterfactruals
[sst_GWp, irfg_GWp, sst_SCCp, irfg_SCCp, ss ] = SolveCounterfactuals(optns,ss,p,zetas) ;

end





