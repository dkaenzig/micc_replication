function irfs = varirf(varEst,invA0,hor,unitShock,cum)
%VARIRF computes all IRFs for the structural VAR(p) with structural impact
% matrix invA0
% Inputs:  -varEst: Estimation results from VAR(p). Relevant component
%                   -B: matrix of VAR coefficients (without constant)
%                   -const: Boolean whether a constant is included
%          -invA0: structural impact matrix (satisfies invA0*invA0' = Sigma
%          -hor: horizon of IRFs 
%          -unitShock: Boolean whether to normalize the IRFs by the SD
%          -cum: 1 x n array of Booleans that indicate for which series to
%                compute cumulative IRFs
% Outputs: -irfs: hor x nvar x nshock matrix of impulse responses
%
% Diego R. Känzig. This version: 17/03/2017

if nargin<5
    cum = false(1,size(varEst.B,1));  % default not cumulative
end
if nargin<4
    unitShock = false;  % default: sd shock
end

% get inputs from varEst
B = varEst.B;
% if varEst.const
B = B(:,1+varEst.nexo:end);
% end

if unitShock
    invA0 = invA0*diag(diag(invA0))^(-1); % normalize A0 to get shocks of 1 (and not a sd)
end

N = size(B,1);
F = varcompan(B);   % get companion form
irfs = zeros(N,N,hor+1);  % preallocate
irfs(:,:,1) = invA0;

% compute IRFs
Fj=F;
for i = 2:hor+1
	irfs(:,:,i) = Fj(1:N,1:N)*invA0;
	Fj=Fj*F;
end

irfs = permute(irfs,[3 1 2]); % permute irfs (shocks in last dimension)

if any(cum)~=0  % compute CIRF if requested
    for i = 1:N
        if cum(i)
            irfs(:,i,:) = cumsum(irfs(:,i,:));
        end
    end
end
    
end

% subfunctions
function F = varcompan(B)
% build matrix F (sparse)
% Output:  F - np x np matrix
% Input: phi - n x np matrix of VAR coeffs (without constant)

[n, pn]=size(B);
F=spalloc(pn,pn,(n+1)*pn);
F(1:n,:)=B;
F(n+1:pn,1:pn-n)=speye(pn-n);
end