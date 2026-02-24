function irf_trans = gettransirf(yirf,xirf,xpath)
%GETTRANSIRF get transitory IRF using Sims' method
% inputs: - yirf: the horizon+1 x 1 IRF of the outcome variable
%         - xirf: the horizon+1 x 1 IRF of the shock variable
%         - xpath: the desired horizon+1 x 1 path for the shock variable
% outputs: - irf_trans: the horizon+1 x 1 IRF of the outcome variable
%                       conditioning on desired shock path
% DK, April, 2024

horizon = size(yirf,1)-1;

% (i) compute the shocks to x that deliver the desired x path 
B = eye(horizon+1); 
for h = 1:horizon
    for i = 1:h
		B(h+1,i) = xirf(h-i+2,1)';
    end
end
epsx = inv(B)*xpath;

% (ii) compute IRF and its covariance matrix wrt the x shocks
shockmat = eye(horizon+1);
for i = 1:horizon+1
    for j = i:horizon+1
		shockmat(i,j) = epsx(j-i+1,1);
    end
end
irf_trans = shockmat'*yirf;

end

