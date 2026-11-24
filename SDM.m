%fn = @fn8_1
%x0 = [ 4 2 -1 ]

function minimo = SDM( fn, x0 )

    iterMax = 10000;
    iterPmax = iterMax;
    iterSmax = iterMax;
    
    tol = 1e-6;
    tolP = tol;
    tolS = tolP/10;
    
    dim = length(x0);
    
    x = x0; % Valor inicial
    iterP = 1;
    [ G, F ] = gradiente( fn, x );
    while( (norm(G)>tolP) && (iterP <= iterPmax) ) 
     
        % NEWTON 1D TRAZA
        alfa = 0; % Valor inicial
        iterS = 1;
        delta = 1;
        while ( (abs(delta)>tolS) && (iterS <= iterSmax) )
            
            AA = DualE(alfa,1,0); % Dual Extendido
            XX = cell(dim,1);
            for i = 1:dim
                XX{i} = x(i) - AA*G(i);
            end
            
            YY = fn( XX );
            
            delta = YY.D/YY.DD;
            
            alfa = alfa - delta;
        
            iterS = iterS + 1;
            
        end
        
        fprintf('iterP = %5d, iterS = %5d, X = [ ', iterP, iterS)
        fprintf('%g ',x)
        fprintf(']\n')
        
        for i = 1:dim
            x(i) = XX{i}.f;
        end
        iterP = iterP+1;
        [ G, F ] = gradiente( fn, x );
        
    end
    
    minimo = x;
        
        
        
end
