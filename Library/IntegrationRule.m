classdef IntegrationRule < handle 
    %UNTITLED Summary of this class goes here
    %   Loool
    
    properties
        gps % array with all the gps 
    end
    
    methods
        
        function obj = setupCubeRule(obj, ngp_xi, ngp_eta, ngp_zeta)
            
            %Get gauess weights and points
            [gpx, gwx] = lgwt(ngp_xi,-1,1);   %gausQuadTable(ngp_xi);
            [gpy, gwy] = lgwt(ngp_eta,-1,1);  %gausQuadTable(ngp_eta);
            [gpz, gwz] = lgwt(ngp_zeta,-1,1); %gausQuadTable(ngp_zeta);
            ngps = ngp_xi*ngp_eta*ngp_zeta;
            
            % allocate array with new Gauss points
            gps(1,ngps) = GaussPoint();
            
            % Set data for each gp
            num = 1;
            for i = 1:ngp_xi
                for j = 1:ngp_eta
                    for k = 1:ngp_zeta
                        local_coords = [gpx(i) gpy(j) gpz(k)];
                        weight = gwx(i) * gwy(j) * gwz(k);
                        %obj.gps{number} = GaussPoint();
                        gps(num).weight = weight;
                        gps(num).number = num;
                        gps(num).local_coords = local_coords;
                        
                        %= GaussPoint(number, local_coords, weight);
                        num = num + 1;
                    end
                end
            end
            
            obj.gps = gps;
            
        end
        
        
    end
    
    end


