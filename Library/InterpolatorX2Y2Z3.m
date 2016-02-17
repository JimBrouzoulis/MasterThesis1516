classdef InterpolatorX2Y2Z3  < Interpolator
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
        
        
    end
    
    methods
        function N = eval_N(obj, localcoords)
            xi = localcoords(1); eta = localcoords(2); zeta = localcoords(3);
            N = [ zeta*(xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2), -zeta*(xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2), -zeta*(xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2), zeta*(xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2), -(xi/2 - 1/2)*(eta/2 - 1/2)*(zeta - 1)*(zeta + 1), (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta - 1)*(zeta + 1), (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta - 1)*(zeta + 1), -(xi/2 + 1/2)*(eta/2 + 1/2)*(zeta - 1)*(zeta + 1), zeta*(xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2), -zeta*(xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2), -zeta*(xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2), zeta*(xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2)];
        end
        
        function [dNdx, detJ] = eval_dNdx(obj, localcoords, ex, ey, ez)
            xi = localcoords(1); eta = localcoords(2); zeta = localcoords(3);

            dNdxi =...
            [(zeta*(eta/2 - 1/2)*(zeta/2 - 1/2))/2,                                        -(zeta*(eta/2 - 1/2)*(zeta/2 - 1/2))/2,                                        -(zeta*(eta/2 + 1/2)*(zeta/2 - 1/2))/2,                                       (zeta*(eta/2 + 1/2)*(zeta/2 - 1/2))/2,                                    -((eta/2 - 1/2)*(zeta - 1)*(zeta + 1))/2,                                   ((eta/2 - 1/2)*(zeta - 1)*(zeta + 1))/2,                                   ((eta/2 + 1/2)*(zeta - 1)*(zeta + 1))/2,                                    -((eta/2 + 1/2)*(zeta - 1)*(zeta + 1))/2,                                       (zeta*(eta/2 - 1/2)*(zeta/2 + 1/2))/2,                                        -(zeta*(eta/2 - 1/2)*(zeta/2 + 1/2))/2,                                        -(zeta*(eta/2 + 1/2)*(zeta/2 + 1/2))/2,                                       (zeta*(eta/2 + 1/2)*(zeta/2 + 1/2))/2;...
              (zeta*(xi/2 - 1/2)*(zeta/2 - 1/2))/2,                                        -(zeta*(xi/2 + 1/2)*(zeta/2 - 1/2))/2,                                        -(zeta*(xi/2 - 1/2)*(zeta/2 - 1/2))/2,                                       (zeta*(xi/2 + 1/2)*(zeta/2 - 1/2))/2,                                    -((xi/2 - 1/2)*(zeta - 1)*(zeta + 1))/2,                                   ((xi/2 + 1/2)*(zeta - 1)*(zeta + 1))/2,                                   ((xi/2 - 1/2)*(zeta - 1)*(zeta + 1))/2,                                    -((xi/2 + 1/2)*(zeta - 1)*(zeta + 1))/2,                                       (zeta*(xi/2 - 1/2)*(zeta/2 + 1/2))/2,                                        -(zeta*(xi/2 + 1/2)*(zeta/2 + 1/2))/2,                                        -(zeta*(xi/2 - 1/2)*(zeta/2 + 1/2))/2,                                       (zeta*(xi/2 + 1/2)*(zeta/2 + 1/2))/2;...
              (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2) + (zeta*(xi/2 - 1/2)*(eta/2 - 1/2))/2, - (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2) - (zeta*(xi/2 + 1/2)*(eta/2 - 1/2))/2, - (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2) - (zeta*(xi/2 - 1/2)*(eta/2 + 1/2))/2, (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2) + (zeta*(xi/2 + 1/2)*(eta/2 + 1/2))/2, - (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta - 1) - (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta + 1), (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta - 1) + (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta + 1), (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta - 1) + (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta + 1), - (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta - 1) - (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta + 1), (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2) + (zeta*(xi/2 - 1/2)*(eta/2 - 1/2))/2, - (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2) - (zeta*(xi/2 + 1/2)*(eta/2 - 1/2))/2, - (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2) - (zeta*(xi/2 - 1/2)*(eta/2 + 1/2))/2, (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2) + (zeta*(xi/2 + 1/2)*(eta/2 + 1/2))/2];    

            %Jacobian transpose
            JT = dNdxi*[ex', ey', ez'];
            detJ=det(JT);
            if (detJ<0)
                error('Jacobian not invertable')
            end

            %Derivatives of x and y
            dNdx = JT\dNdxi;
        end
        
    end
    
end

