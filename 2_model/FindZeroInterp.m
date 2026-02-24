function [ Value , iL , iU , frac ] = FindZeroInterp( Function , Grid )

% finds the value, on the grid Grid, that sets the function Function to 0

% Arguments must be column vectors
if size(Function,1) == 1
    Function = Function' ;
end

if size(Grid,1) == 1
    Grid = Grid' ;
end

N = size(Function,1) ;

% Find closest point to zero
[~,imin] = min( Function.^2 ) ;

% Determine whether function is locally increasing or decreasing
imin = min( max( imin , 2 ) , N - 1 ) ;

if Function(imin+1) > Function(imin-1) % Function locally increasing
    
    if Function(imin) > 0
        imin = imin  - 1 ;
    end
    iL = imin; iU = imin + 1 ;
    frac = ( 0 - Function(iL) ) / ( Function(iU) - Function(iL) ) ;
    Value = ( 1 - frac ) * Grid(iL) + frac * Grid(iU) ;
  
elseif Function(imin+1) < Function(imin-1) % Function locally decreasing
    
    if Function(imin) < 0
        imin = imin  - 1 ;
    end
    iL = imin; iU = imin + 1 ;
%     frac = Function(iU) / ( Function(iU) - Function(iL) ) ;
    frac = ( 0 - Function(iL) ) / ( Function(iU) - Function(iL) ) ;
    Value = ( 1 - frac ) * Grid(iL) + frac * Grid(iU) ;
   
else 
    % disp('WARNING: Function locally flat, using default index for interpolation')
    iL = imin + 1 ; iU = imin + 2 ;
    frac = 0.5 ;
    Value = ( 1 - frac ) * Grid(iL) + frac * Grid(iU) ;

end

end

