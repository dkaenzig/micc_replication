function bootData = varbootstrap(varEst,initVal,dataExo,nsim,boottype)
%VARBOOTSTRAP generates nsim artificial time series by drawing with replacement
% from the shocks and feeding them to the estimated VAR(p)
%
% Inputs: -varEst: estimation results of VAR(p). Relevant components:
%                  -B: estimated VAR coefficients
%                  -const: Boolean whether to include a constant
%                  -U: residuals of the VAR
%         -initVal: p x nvar matrix of initial values
%         -nsim: number of simulations
%         -boottype: String specifying how to bootstrap: standard bootstrap
%                    by drawing from residuals with replacemend: 'bs', wild bootstrap
%                    based on Rademacher distribution: 'wild'
% Output: -bootData: bootstrapped data
%
% Diego R. Känzig. This version: 30/03/2017

if nargin<5
    boottype = 'bs';
end

T = varEst.T;
n = varEst.n;
p = varEst.p;

bootData = zeros(T+p,n,nsim); % 3rd index: which boostrap sample

for j=1:nsim
    
    if strcmp(boottype,'bs')
        index = randsample(1:T,T,true)';    % test code using (1:T)';
        Uboot = varEst.U(index,:);
    elseif strcmp(boottype,'wild')          % use wild bootstrap based on Rademacher distribution
        rr = 1-2*(rand(T,1)>0.5);
        Uboot = (varEst.U).*(rr*ones(1,n));
    end
    
    bootData(1:p,:,j) = initVal; %initial values of y, same for all j
    for i = p+1:T+p
        bootData(i,:,j)= varEst.B*[dataExo(i,:)'; vec(fliplr(bootData(i-p:i-1,:,j)'))] ...
                         + Uboot(i-p,:)'; % bootstrap
    end

end

end
