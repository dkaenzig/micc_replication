

%%% parameters for SCC calcualtion

% world GDP and consumption in 2024 international dollars, from world Bank https://data.worldbank.org/indicator/NY.GDP.MKTP.PP.CD
p.WorldGDP  = 198.7 * 1e12 ; 
p.WorldCons = 0.8 * p.WorldGDP ; % assumes 20% goes in investment

% size of the carbon pulse corresponding to temperature path below
p.Gt100    = 100 * 1e9 * 3.667 ;
    % 3.667 adjusts reporting in GtC vs. GtCO2 in Joos et al. (2013) / Dietz et al. (2021)

% scale of carbon pulse in fraction of 100 GtCO2
% adjust downward to consider truly marginal change and not pick up model concavity
p.CO2size.current = 1 ; % default, will be set to one of the two values below so that temperature IRF is of comparable magnitude
p.CO2size.global  = 0.1 ; 
p.CO2size.local   = 8 * p.CO2size.global ; % to have comparable implied GDP path for global and local temperature damages


% temperature path from Dietz et al. (2021)
% response in Celsius to a 100 GtC pulse
% fit parametric form: T = a*( exp(-b*t) - exp(-c*t) ) + d*( 1 - exp(-f*t) ) ;

% adjustable scale of temperature pulse
p.adj    = 1 ; 

% temperature values at time horizons
y0   = 0 ;
y5   = 0.18 ;
y10  = 0.215 ;
y20  = 0.2 ;
y40  = 0.18 ;
y200 = 0.165 ;

ys = [y5 y10 y20 y40 y200 ]' ;
ts = [ 5 10  20  40  200 ]' ;

% fit
g   = @(x) x(1) * ( exp( -x(2) * ts ) - exp( -x(3) * ts ) ) + x(4)*( 1 - exp( -x(5) * ts ) ) ; % temperature path
obj = @(x) 1e4 * sum( ( g(x) - ys ).^2 ./ ts ) ; % weight by inverse time to capture initial behavior
x0  = [ 0.1 0.1 0.1 0.1 y200 ]' ; % initial guess

% run optimizer
options = optimoptions('fmincon', 'Display', 'off') ;
xd = fmincon(obj , ... objective function
             x0 , ... starting point
             [],[],[],[], ... generaly inequality and equality constraints, not used here
             [ 0 ; 0 ; 0 ; 0 ; 0 ] , [ 1 ; 10 ; 10 ; 1 ; 10 ] ,...  lower and upper bound on parameters
             [] , options ) ; % nonlinear constraints and options
p.xD = xd ; % save estimate





