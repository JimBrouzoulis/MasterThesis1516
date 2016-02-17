classdef InterpolatorX2Y2Z2 < Interpolator
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function N = eval_N(obj, localcoords)
            xi = localcoords(1); eta = localcoords(2); zeta = localcoords(3);
            N = [ -(xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2), (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2), (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2), -(xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2), (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2), -(xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2), -(xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2), (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2)];
        end
        
        function [dNdx, detJ] = eval_dNdx(obj, lcoords, ex, ey, ez)
            xi = lcoords(1); eta = lcoords(2); zeta = lcoords(3);
            dNdxi =...
            [ -((eta/2 - 1/2)*(zeta/2 - 1/2))/2, ((eta/2 - 1/2)*(zeta/2 - 1/2))/2, ((eta/2 + 1/2)*(zeta/2 - 1/2))/2, -((eta/2 + 1/2)*(zeta/2 - 1/2))/2, ((eta/2 - 1/2)*(zeta/2 + 1/2))/2, -((eta/2 - 1/2)*(zeta/2 + 1/2))/2, -((eta/2 + 1/2)*(zeta/2 + 1/2))/2, ((eta/2 + 1/2)*(zeta/2 + 1/2))/2;...
             -((xi/2 - 1/2)*(zeta/2 - 1/2))/2, ((xi/2 + 1/2)*(zeta/2 - 1/2))/2, ((xi/2 - 1/2)*(zeta/2 - 1/2))/2, -((xi/2 + 1/2)*(zeta/2 - 1/2))/2, ((xi/2 - 1/2)*(zeta/2 + 1/2))/2, -((xi/2 + 1/2)*(zeta/2 + 1/2))/2, -((xi/2 - 1/2)*(zeta/2 + 1/2))/2, ((xi/2 + 1/2)*(zeta/2 + 1/2))/2;...
             -((xi/2 - 1/2)*(eta/2 - 1/2))/2, ((xi/2 + 1/2)*(eta/2 - 1/2))/2, ((xi/2 - 1/2)*(eta/2 + 1/2))/2, -((xi/2 + 1/2)*(eta/2 + 1/2))/2, ((xi/2 - 1/2)*(eta/2 - 1/2))/2, -((xi/2 + 1/2)*(eta/2 - 1/2))/2, -((xi/2 - 1/2)*(eta/2 + 1/2))/2, ((xi/2 + 1/2)*(eta/2 + 1/2))/2];
            
            %Jacobian transpose (but ignore middle point)
            JT=dNdxi*[ex', ey', ez'];
            detJ=det(JT);
            if (detJ<0)
                error('Jacobian not invertable')
            end

            %Derivatives of x and y
            dNdx = JT\dNdxi;
        end
     
    end
    
end

