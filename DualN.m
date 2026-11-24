classdef DualN
   
    properties
        
        f % Valor de la Funcion
        D % Valor de la Derivada de la Funcion
        
    end
    
    methods
        
        % Constructor
        function h = DualN(a,b)
            
            if nargin == 0
                h.f = 0;
                h.D = 0;
            elseif nargin == 1 % Constante, derivada cero
                h.f = a;
                h.D = 0;
            else
                h.f = a;
                h.D = b;
            end
            
        end
        
        % Cambio de signo (multiplicacion por -1)
        function h = uminus(u)
            f = -u.f;
            D = -u.D;
            h = DualN(f,D);
        end
        
        % No Cambio de signo (multiplicacion por +1)
        function h = uplus(u)
            h = u;
        end
        
        % Suma de Derivadas
        function h = plus(u,v)
            if ~isa(u,'DualN') %u es scalar, v es Dual 
                f = u + v.f;
                D = v.D; % Regla de Derivacio 
            elseif ~isa(v,'DualN') %v es scalar, u es Dual
                f = u.f + v;
                D = u.D;  
            else % ambos, u y v son Duales
                f = u.f + v.f;
                D = u.D + v.D; % Regla de Derivacion
            end
            h = DualN( f, D );
        end
        
        % Resta de Derivadas
        function h = minus(u,v)
            h = u + (-v);
        end
        
        % Producto de Derivadas
        function h = mtimes(u,v)
            if ~isa(u,'DualN') %u es scalar, v es Dual 
                f = u * v.f;
                D = u * v.D; % Regla de Derivacion
            elseif ~isa(v,'DualN') %v es scalar, u es Dual
                f = u.f * v;
                D = v * u.D  ; % Regla de Derivacion  
            else % ambos, u y v son Duales
                f = u.f * v.f;
                D = u.f * v.D + v.f * u.D  ; % Regla de Derivacion
            end
            h = DualN( f, D );
        end
        
        % Exponencial 
        function h = exp( u )
            f = exp( u.f );
            D = exp( u.f )*u.D;
            h = DualN( f, D );
        end
        
        % Logaritmo  
        function h = log( u )
            f = log( u.f );
            D = u.D / u.f;
            h = DualN( f, D );
        end
        
        % Potencia
        function h = mpower(u,v)
            if ~isa(u,'DualN') %u es scalar, v es Dual 
                f = u ^ v.f;
                D = ( u ^ v.f ) * log(u) * v.D;
                h = DualN( f, D );
            elseif ~isa(v,'DualN') %v es scalar, u es Dual; x^n
                f = u.f ^ v;
                D = v * ( u.f^(v-1) ) * u.D;
                h = DualN( f, D );
            else % ambos, u y v son Duales
                h = exp( v * log(u) ); 
            end
            
        end
        
        % Seno  
        function h = sin( u )
            f = sin( u.f );
            D = cos( u.f )*u.D;
            h = DualN( f, D );
        end
        
        % Coseno  
        function h = cos( u )
            f = cos( u.f );
            D = -sin( u.f )*u.D;
            h = DualN( f, D );
        end
        
        % Raiz Cuadrada  
        function h = sqrt( u )
            f = sqrt( u.f );
            D = u.D/(2*sqrt( u.f ));
            h = DualN( f, D );
        end
        
        
        
    end
    
    
end