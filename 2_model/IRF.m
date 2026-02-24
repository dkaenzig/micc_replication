function irf = IRF(p,n,ss,irf,tra)

%%% computes transitional dynamics

K     = zeros(p.Nt,1) ;
E     = zeros(p.Nt,1) ;
I     = zeros(p.Nt,1) ;
V     = zeros(p.Nt,1) ;
W     = zeros(p.Nt,1) ;

K(1) = irf.K0 ;

for t=1:p.Nt

    [e,sav,v] = evalPolicy(n, K(t), tra.E(:,t), tra.sav(:,t), tra.V(:,t) ) ;
    E(t)    = e ;
    V(t)    = v ;
    I(t)    = sav + irf.Delta(t) * K(t) ;
    
    if p.gamma ~= 1
        W(t) = ( v/ss.Viss )^(1/(1-p.gamma))- 1 ;
    else
        W(t) = exp( p.rho*(v-ss.Viss) ) - 1 ;
    end

    if t < p.Nt
        K(t+1)  = K(t) + sav * p.dt ; 
    end

end

% output in percent deviations or pp
irf.V       = V ;
irf.Kp      = K/ss.K-1 ;
irf.Ep      = E/ss.E-1 ;
irf.Zp      = irf.zhat ;
irf.Delta   = irf.Delta ;
irf.Yp      = exp(irf.zhat ) .* ( K/ss.K ).^p.alpha - 1 ;
irf.Ip      = I / ss.I -1 ;
irf.welfare = W ;

end