function L = ConvolutionMatrix(p,T)
    
    %%% matrix to convert from transitory to persistent; corresponds to convolution operation
    
    L  = zeros(p.Nt) ;
    for i=1:p.Nt
        L(i,i) = 1 ;
        if i > 1
            for j=1:(i-1)
                L(i,j) = T(i-j+1) ;
            end
        end
    
    end

end