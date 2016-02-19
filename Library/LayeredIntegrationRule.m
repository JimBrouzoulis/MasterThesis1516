classdef LayeredIntegrationRule < handle
    %UNTITLED Summary of this class goes here
    %   Loool
    
    properties
        %Array off all IR for all plies
        irs = IntegrationRule;
    end
    
    methods
        
        function obj = LayeredIntegrationRule(nlay, ngp_xi, ngp_eta, ngp_zeta)
            
            %Create one IR for each layer
            for il = 1:nlay
                obj.irs(il) = IntegrationRule;
                obj.irs(il).setupCubeRule(ngp_xi, ngp_eta ,ngp_zeta);
            end
            
        end
        
        function egp = getElementGaussCoordinate(obj, ilay, l_z_coord, layZCoords)         
            zc = l_z_coord; %zc = z-coordinate 

            %Shape function in 1D
            N1D = [(-zc + 1)/2, (zc + 1)/2];

            %interpolate to get element gauss point coord
            egp = N1D(1)*layZCoords(1) + N1D(2)*layZCoords(2);
        end
    end
    
end


