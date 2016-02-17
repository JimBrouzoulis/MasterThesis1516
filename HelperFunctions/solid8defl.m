function [ disp_x,disp_y,disp_z ] = deflectionField( xsi,eta,zeta,Ex,Ey,Ez,Ed)

% Ex =    [x1, x1, ...
%          x2, x2, ...
%           :   :  ...
%          x8, x8  ... ];  [8 x nel]

nx = length(xsi);
ny = length(eta);
nz = length(zeta);
nel = size(Ex,2);

disp_x = zeros(nx,ny,nz);
disp_y = zeros(nx,ny,nz);
disp_z = zeros(nx,ny,nz);

for xx = 1:nx
    for yy = 1:ny
        for zz = 1:nz
           
            for inel = 1:nel
               
                cEx = Ex(:,inel);
                cEy = Ey(:,inel);
                cEz = Ez(:,inel);
                shp = alphaShape(cEx,cEy,cEz);
                shp.Alpha = shp.Alpha*10;
                isInShape = inShape(shp,xsi(xx),eta(yy),zeta(zz));
                if(isInShape)
                    meanX = mean(cEx);
                    meanY = mean(cEy);
                    meanZ = mean(cEz);
                    Nxieta = shapequad(xsi(xx)-meanX,eta(yy)-meanY, zeta(zz)-meanZ);
                    Nxy = [Nxieta(1)*eye(3),Nxieta(2)*eye(3),Nxieta(3)*eye(3) ,Nxieta(4)*eye(3),...
                           Nxieta(5)*eye(3) ,Nxieta(6)*eye(3), Nxieta(7)*eye(3) , Nxieta(8)*eye(3)];
                    
                    ae = Ed(:,inel);
                    u_xyz = Nxy*ae;
                    
                    disp_x(xx,yy,zz) = u_xyz(1);
                    disp_y(xx,yy,zz) = u_xyz(2);
                    disp_z(xx,yy,zz) = u_xyz(3);
                    
                    break;
                end
            end
            if(inel == nel)
               disp('point not find in any elemtn'); 
            end
        end
    end
end

end


function [N,B]=shapequad(xi,eta,zeta)

N(1) = 0.125*(1 - xi)*(1 - eta)*(1 - zeta);
N(2) = 0.125*(1 + xi)*(1 - eta)*(1 - zeta);
N(3) = 0.125*(1 + xi)*(1 + eta)*(1 - zeta);
N(4) = 0.125*(1 - xi)*(1 + eta)*(1 - zeta);
N(5) = 0.125*(1 - xi)*(1 - eta)*(1 + zeta);
N(6) = 0.125*(1 + xi)*(1 - eta)*(1 + zeta);
N(7) = 0.125*(1 + xi)*(1 + eta)*(1 + zeta);
N(8) = 0.125*(1 - xi)*(1 + eta)*(1 + zeta);

B = [ -((eta - 1)*(zeta - 1))/8, ((eta - 1)*(zeta - 1))/8, -((eta + 1)*(zeta - 1))/8, ((eta + 1)*(zeta - 1))/8, ((eta - 1)*(zeta + 1))/8, -((eta - 1)*(zeta + 1))/8, ((eta + 1)*(zeta + 1))/8, -((eta + 1)*(zeta + 1))/8;...
    -(xi/8 - 1/8)*(zeta - 1),  (xi/8 + 1/8)*(zeta - 1),  -(xi/8 + 1/8)*(zeta - 1),  (xi/8 - 1/8)*(zeta - 1),  (xi/8 - 1/8)*(zeta + 1),  -(xi/8 + 1/8)*(zeta + 1),  (xi/8 + 1/8)*(zeta + 1),  -(xi/8 - 1/8)*(zeta + 1);...
    -(xi/8 - 1/8)*(eta - 1),   (xi/8 + 1/8)*(eta - 1),   -(xi/8 + 1/8)*(eta + 1),   (xi/8 - 1/8)*(eta + 1),   (xi/8 - 1/8)*(eta - 1),   -(xi/8 + 1/8)*(eta - 1),   (xi/8 + 1/8)*(eta + 1),   -(xi/8 - 1/8)*(eta + 1)];
end

