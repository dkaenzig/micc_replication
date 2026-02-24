
%%% externally set model parameters
% parameters at annual frequency

% househod
p.gamma = 1 ;           % risk aversion.

% utility function
if p.gamma ~= 1
        u = @(E) ( E.^( 1 - p.gamma ) ) / ( 1 - p.gamma ) ; % exclude the -1 for welfare gains formula
else
        u = @(E) log(E) ;
end
up  = @(E) E.^(-p.gamma) ;
up1 = @(E) E.^(-1/p.gamma) ;

% production
p.Z      = 1 ;
p.alpha  = 0.33 ;
p.Delta  = 0.08 ;



