classdef DualH
   
    properties
        
        f % Valor de la Funcion
        D1 % Valor de la Derivada con respecto a X1 de la Funcion
        D2 % Valor de la Derivada con respecto a X2 de la Funcion
        DD % Valor de la Segunda Derivada Cruzada de la Funcion
        
    end % properties
   
    methods
        
        % Constructor
        function h = DualH(a,b,c,d)
            
            tol_cero = 2*eps('double');
            
            if nargin == 0
                h.f = 0;
                h.D1 = 0;
                h.D2 = 0;
                h.DD = 0;
            elseif nargin == 1 % Constante, derivada cero
                h.f = a;
                h.D1 = 0;
                h.D2 = 0;
                h.DD = 0;
            elseif nargin == 3 
                %condicion
                h.f = a;
                h.D1 = b; % verificar 1
                h.D2 = c; % verificar 1
                h.DD = 0;
            else
                h.f = a;
                h.D1 = b;
                h.D2 = c;
                h.DD = d;
            end
            
        end
        
        % No Cambio de signo (multiplicacion por +1)
        function h = uplus(u)
            h = u;
        end
            
        % Cambio de signo (multiplicacion por -1)
        function h = uminus(u)
            f = -u.f;
            D1 = -u.D1;
            D2 = -u.D2;
            DD = -u.DD;
            h = DualH(f,D1,D2,DD);
        end
        
        % Suma de Derivadas
        function h = plus(u,v)
            if ~isa(u,'DualH') %u es scalar, v es Dual 
                f = u + v.f;
                D1 = v.D1; % Regla 1° Derivacion con respecto a X1
                D2 = v.D2; % Regla 1° Derivacion con respecto a X2
                DD = v.DD; % Regla 2° Derivacion Cruzada
            elseif ~isa(v,'DualH') %v es scalar, u es Dual
                f = u.f + v;
                D1 = u.D1; % Regla 1° Derivacion con respecto a X1
                D2 = u.D2; % Regla 1° Derivacion con respecto a X2
                DD = u.DD; % Regla 2° Derivacion Cruzada  
            else % ambos, u y v son Duales
                f = u.f + v.f;
                D1 = u.D1 + v.D1; % Regla 1° Derivacion con respecto a X1
                D2 = u.D2 + v.D2; % Regla 1° Derivacion con respecto a X2
                DD = u.DD + v.DD; % Regla 2° Derivacion Cruzada
            end
            h = DualH( f, D1, D2, DD );
        end
        
        % Resta de Derivadas
        function h = minus(u,v)
            h = u + (-v);
        end
        
        % Producto de Derivadas
        function h = mtimes(u,v)
            if ~isa(u,'DualH') %u es scalar, v es Dual 
                f = u * v.f;
                D1 = u * v.D1;  % Regla 1° Derivacion con respecto a X1
                D2 = u * v.D2;  % Regla 1° Derivacion con respecto a X2
                DD = u * v.DD;  % Regla 2° Derivacion Cruzada
            elseif ~isa(v,'DualH') %v es scalar, u es Dual
                f = v * u.f;
                D1 = v * u.D1;  % Regla 1° Derivacion con respecto a X1
                D2 = v * u.D2;  % Regla 1° Derivacion con respecto a X2
                DD = v * u.DD;  % Regla 2° Derivacion Cruzada
            else % ambos, u y v son Duales
                f = u.f * v.f;
                D1 = v.f * u.D1 + u.f * v.D1;  % Regla 1° Derivacion con respecto a X1
                D2 = v.f * u.D2 + u.f * v.D2;  % Regla 1° Derivacion con respecto a X2
                DD = u.D1 * v.D2 + v.D1 * u.D2 + v.f * u.DD + u.f * v.DD;  % Regla 2° Derivacion Cruzada
            end
            h = DualH( f, D1, D2, DD );
        end
        
        % Producto de Potencias
        function h = mpower(u,v)
            if ~isa(u,'DualH') %u es scalar, v es Dual 
                error()
            elseif ~isa(v,'DualH') %v es scalar, u es Dual
                f = u.f ^ v;
                D1 = v * u.f^(v-1) * u.D1;
                D2 = v * u.f^(v-1) * u.D2;
                DD = v * ( (v-1)*u.f^(v-2) * u.D1 * u.D2 + u.f^(v-1) * u.DD );
            else % ambos, u y v son Duales
                error()
            end
            h = DualH( f, D1, D2, DD );
        end
        
        
        
    end % methods
    
end % classdef