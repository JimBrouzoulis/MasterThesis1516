classdef ElementProperties < handle
    %ElementProperties Storage for both material data and cross sectional
    %data (could be separated into two classes if one wants to)
    
    
    properties
        
        thickness
        angles
        int_coordsG % Coordinates for interfaces between lamina
        int_coordsL
        D % Material stiffness
        lamZCoordsL
        lamZCoordsG
        nLam
        Dmatrices
    end
    
    methods
        function obj = setup(obj, thickness, angles, D)
            obj.thickness = thickness;
            obj.angles = pi/180*angles;  
            hlam = obj.thickness/length(angles); % Thickness of a lamina
            
            % TODO: should this be local coords or in global? should then be
            % set with respect to the global z-coords
            obj.int_coordsG = 0:hlam:obj.thickness;
            
            
            %Variable lamCoords are in the global system, convert to the
            %element local system, from [-1 to 1]
            ic = obj.int_coordsG;
            obj.int_coordsL = ...
                (2*obj.int_coordsG - (ic(end) + ic(1)) ) / (ic(end)-ic(1));
            obj.nLam = length(angles);
            
           
            %Define D-matrices
            for id = 1:obj.nLam
                T1 = rotmat(obj.angles(id),1); 
                T2 = rotmat(obj.angles(id),2);
                obj.Dmatrices(:,:,id) = T1^-1 * D * T2;
            end
            
            
        end
        
        
    end
    
end

