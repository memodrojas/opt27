classdef DualE
   
    properties
        
        f % Valor de la Funcion
        D % Valor de la Derivada de la Funcion
        DD % Valor de la Segunda Derivada de la Funcion
        
    end
   
    methods
        
        % Constructor
        function h = DualE(a,b,c)
            
            tol_cero = 2*eps('double');
            
            if nargin == 0
                h.f = 0;
                h.D = 0;
                h.DD = 0;
            elseif nargin == 1 % Constante, derivada cero
                h.f = a;
                h.D = 0;
                h.DD = 0;
            elseif nargin == 2 % Variable independiente (derivada 1)
                if abs( b-1 ) < tol_cero
                    h.f = a;
                    h.D = 1;
                    h.DD = 0;
                else
                    error('DualE: 1° derivada no es 1.')
                end
            else
                h.f = a;
                h.D = b;
                h.DD = c;
            end
            
        end
       
        % Cambio de signo (multiplicacion por -1)
        function h = uminus(u)
            f = -u.f;
            D = -u.D;
            DD = -u.DD;
            h = DualE(f,D,DD);
        end
        
        % No Cambio de signo (multiplicacion por +1)
        function h = uplus(u)
            h = u;
        end
        
        % Suma de Derivadas
        function h = plus(u,v)
            if ~isa(u,'DualE') %u es scalar, v es Dual 
                f = u + v.f;
                D = v.D; % Regla 1° Derivacion
                DD = v.DD; % Regla 2° Derivacion
            elseif ~isa(v,'DualE') %v es scalar, u es Dual
                f = u.f + v;
                D = u.D; % Regla 1° Derivacion
                DD = u.DD; % Regla 2° Derivacion  
            else % ambos, u y v son Duales
                f = u.f + v.f;
                D = u.D + v.D; % Regla 1° Derivacion
                DD = u.DD + v.DD; % Regla 2° Derivacion
            end
            h = DualE( f, D, DD );
        end
        
        % Resta de Derivadas
        function h = minus(u,v)
            h = u + (-v);
        end
        
        % Producto de Derivadas
        function h = mtimes(u,v)
            if ~isa(u,'DualE') %u es scalar, v es Dual 
                f = u * v.f;
                D = u * v.D; % Regla 1° Derivacion
                DD = u * v.DD;  % Regla 2° Derivacion
            elseif ~isa(v,'DualE') %v es scalar, u es Dual
                f = u.f * v;
                D = v * u.D  ; % Regla 1° Derivacion  
                DD = v * u.DD;  % Regla 2° Derivacion
            else % ambos, u y v son Duales
                f = u.f * v.f;
                D = u.f * v.D + v.f * u.D  ; % Regla 1° Derivacion
                DD = u.f * v.DD + u.D * v.D + v.f * u.DD + v.D * u.D;  % Regla 2° Derivacion
            end
            h = DualE( f, D, DD );
        end
        
        % Exponencial 
        function h = exp( u )
            f = exp( u.f );
            D = exp( u.f )*u.D;
            DD = exp( u.f )*u.DD + exp( u.f )*u.D*u.D
            h = DualE( f, D, DD );
        end
        
        % Logaritmo  
        function h = log( u )
            f = log( u.f );
            D = u.D / u.f;
            DD = ( u.f * u.DD - u.D * u.D ) / u.f^2; 
            h = DualE( f, D, DD );
        end
               
        % Potencia
        function h = mpower(u,v)
            if ~isa(u,'DualE') %u es scalar, v es Dual 
                f = u ^ v.f;
                D = ( u ^ v.f ) * log(u) * v.D;
                DD = log(u) * ( ( u ^ v.f ) * log(u) * v.D * v.D + ( u ^ v.f ) * v.DD );
                h = DualE( f, D, DD );
            elseif ~isa(v,'DualE') %v es scalar, u es Dual; x^n
                f = u.f ^ v;
                D = v * ( u.f^(v-1) ) * u.D;
                DD = v * ( (v-1)*u.f^(v-2)*u.D * u.D + ( u.f^(v-1) ) * u.DD  );
                h = DualE( f, D, DD );
            else % ambos, u y v son Duales
                h = exp( v * log(u) ); 
            end
            
        end
        
        % Seno  
        function h = sin( u )
            f = sin( u.f );
            D = cos( u.f )*u.D;
            DD = -sin( u.f )*u.D * u.D + cos( u.f ) * u.DD; 
            h = DualE( f, D, DD );
        end
        
        % Coseno  
        function h = cos( u )
            f = cos( u.f );
            D = -sin( u.f )*u.D;
            DD = - sin( u.f )*u.DD - cos( u.f )*u.D^2;
            h = DualE( f, D, DD );
        end
        
        % Raiz Cuadrada  
        function h = sqrt( u )
            f = sqrt( u.f );
            D = u.D/(2*sqrt( u.f ));
            DD = ( 2*sqrt( u.f )*u.DD - u.D*(2*u.D/(2*sqrt(u.f))) )/(2*sqrt( u.f ))^2;
            h = DualE( f, D, DD );
        end
        
    end
    
end