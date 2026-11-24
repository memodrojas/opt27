function raiz = NewtonR1D( fn, x0 )

    tol = 1e-6;
    MAXITER = 10000;
    
    x = x0;

    iter = 1;
    delta = 1;
    while ( abs(delta) > tol ) && ( iter <= MAXITER )
    
        XX = DualN( x, 1 );
        YY = fn( XX );
        
        delta = YY.f / YY.D;
        
        x = x - delta;
        
        fprintf( 'iter = %d, delta = %g\n ', iter, delta )
        iter = iter + 1;
    end
    
    raiz = x;