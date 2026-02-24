function [E,Sav,V] = evalPolicy(n,K,e,sav,v)

% evaluates consumption, savings and value functions
[ ~ , il , iu , f ] = ...
        FindZeroInterp( K - n.K , n.K ) ;

E   = ( 1 - f ) * e(il)   + f * e(iu) ;
Sav = ( 1 - f ) * sav(il) + f * sav(iu) ;
V   = ( 1 - f ) * v(il)   + f * v(iu) ;

end