function laggedvar = panellagmatrix(var,panelvar,lags)
%PANELLAGMATRIX construct lag matrix for panel variables
% inputs: - var: NT x K matrix with panel data
%         - panelvar: panel identifier
%         - lags: vector with lags
% outputs: - laggedvar: NT x K*lags matrix with lagged variables
% DK, April, 2024

panelvars = unique(panelvar);
N = size(panelvars,1);

% preallocate
laggedvar = nan(size(var,1),size(var,2)*length(lags));

% create lags by panel unit
for ii = 1:N
    pvar = panelvars(ii);
    laggedvar(strcmp(panelvar,pvar),:) = lagmatrix(var(strcmp(panelvar,pvar),:),lags);
end

