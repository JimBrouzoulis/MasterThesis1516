classdef Interpolator 
    %INTERPOLATOR Base class
    %   Detailed explanation goes here
    
    properties
    end
    
    methods
        
        function N = eval_N(obj, lcoords)
            error('This is a abstract function')
        end
        
        function [dNdx, detJ] = eval_dNdx(obj, lcoords, ex, ey, ez)
            error('This is a abstract function')
        end
        
        function [gc] = eval_globalCoords(obj, localCoords,ex,ey,ez)
            N = obj.eval_N(localCoords);
            gc = N*[ex',ey',ez'];
        end
        
        function N = createNmatrix(obj, N_vec, dim)
            N = zeros(dim, length(N_vec));
            
            for i = 1:length(N_vec)
                for j = 1:dim
                    N(j, (i-1)*dim + j ) = N_vec(i);
                end
            end
            
        end
        
        function dofs = getDofsFromNodes(obj, component, nodes, dim)
           
            dofs = nodes*dim - (dim-component);
            
        end
        
    end
    
end

