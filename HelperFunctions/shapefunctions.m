function [ N,B ] = shapefunctions( xi,eta,zeta, opt)
%Options: 1 = linear, 2 = quadratic, 3 = ... and so on... 
%(See explanations below)
%Add according to your need

if(nargin <4)
    opt = 1;
end

%Linear at corder nodes of cube
if opt == 1
% N(1) = 0.125*(1 - xi)*(1 - eta)*(1 - zeta);
% N(2) = 0.125*(1 + xi)*(1 - eta)*(1 - zeta);
% N(3) = 0.125*(1 + xi)*(1 + eta)*(1 - zeta);
% N(4) = 0.125*(1 - xi)*(1 + eta)*(1 - zeta);
% N(5) = 0.125*(1 - xi)*(1 - eta)*(1 + zeta);
% N(6) = 0.125*(1 + xi)*(1 - eta)*(1 + zeta);
% N(7) = 0.125*(1 + xi)*(1 + eta)*(1 + zeta);
% N(8) = 0.125*(1 - xi)*(1 + eta)*(1 + zeta);
% 
% B = [ -((eta - 1)*(zeta - 1))/8, ((eta - 1)*(zeta - 1))/8, -((eta + 1)*(zeta - 1))/8, ((eta + 1)*(zeta - 1))/8, ((eta - 1)*(zeta + 1))/8, -((eta - 1)*(zeta + 1))/8, ((eta + 1)*(zeta + 1))/8, -((eta + 1)*(zeta + 1))/8;...
%     -(xi/8 - 1/8)*(zeta - 1),  (xi/8 + 1/8)*(zeta - 1),  -(xi/8 + 1/8)*(zeta - 1),  (xi/8 - 1/8)*(zeta - 1),  (xi/8 - 1/8)*(zeta + 1),  -(xi/8 + 1/8)*(zeta + 1),  (xi/8 + 1/8)*(zeta + 1),  -(xi/8 - 1/8)*(zeta + 1);...
%     -(xi/8 - 1/8)*(eta - 1),   (xi/8 + 1/8)*(eta - 1),   -(xi/8 + 1/8)*(eta + 1),   (xi/8 - 1/8)*(eta + 1),   (xi/8 - 1/8)*(eta - 1),   -(xi/8 + 1/8)*(eta - 1),   (xi/8 + 1/8)*(eta + 1),   -(xi/8 - 1/8)*(eta + 1)];

N = [ -(xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2), (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2), (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2), -(xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2), (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2), -(xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2), -(xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2), (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2)];
 
B =[ -((eta/2 - 1/2)*(zeta/2 - 1/2))/2, ((eta/2 - 1/2)*(zeta/2 - 1/2))/2, ((eta/2 + 1/2)*(zeta/2 - 1/2))/2, -((eta/2 + 1/2)*(zeta/2 - 1/2))/2, ((eta/2 - 1/2)*(zeta/2 + 1/2))/2, -((eta/2 - 1/2)*(zeta/2 + 1/2))/2, -((eta/2 + 1/2)*(zeta/2 + 1/2))/2, ((eta/2 + 1/2)*(zeta/2 + 1/2))/2;...
 -((xi/2 - 1/2)*(zeta/2 - 1/2))/2, ((xi/2 + 1/2)*(zeta/2 - 1/2))/2, ((xi/2 - 1/2)*(zeta/2 - 1/2))/2, -((xi/2 + 1/2)*(zeta/2 - 1/2))/2, ((xi/2 - 1/2)*(zeta/2 + 1/2))/2, -((xi/2 + 1/2)*(zeta/2 + 1/2))/2, -((xi/2 - 1/2)*(zeta/2 + 1/2))/2, ((xi/2 + 1/2)*(zeta/2 + 1/2))/2;...
 -((xi/2 - 1/2)*(eta/2 - 1/2))/2, ((xi/2 + 1/2)*(eta/2 - 1/2))/2, ((xi/2 - 1/2)*(eta/2 + 1/2))/2, -((xi/2 + 1/2)*(eta/2 + 1/2))/2, ((xi/2 - 1/2)*(eta/2 - 1/2))/2, -((xi/2 + 1/2)*(eta/2 - 1/2))/2, -((xi/2 - 1/2)*(eta/2 + 1/2))/2, ((xi/2 + 1/2)*(eta/2 + 1/2))/2];
 



%Quadratic through thickness, linear in plane. 12 nodes total
elseif opt == 2
    
   N = [ zeta*(xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2), -zeta*(xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2), -zeta*(xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2), zeta*(xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2), -(xi/2 - 1/2)*(eta/2 - 1/2)*(zeta - 1)*(zeta + 1), (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta - 1)*(zeta + 1), (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta - 1)*(zeta + 1), -(xi/2 + 1/2)*(eta/2 + 1/2)*(zeta - 1)*(zeta + 1), zeta*(xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2), -zeta*(xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2), -zeta*(xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2), zeta*(xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2)];
   B =[                                       (zeta*(eta/2 - 1/2)*(zeta/2 - 1/2))/2,                                        -(zeta*(eta/2 - 1/2)*(zeta/2 - 1/2))/2,                                        -(zeta*(eta/2 + 1/2)*(zeta/2 - 1/2))/2,                                       (zeta*(eta/2 + 1/2)*(zeta/2 - 1/2))/2,                                    -((eta/2 - 1/2)*(zeta - 1)*(zeta + 1))/2,                                   ((eta/2 - 1/2)*(zeta - 1)*(zeta + 1))/2,                                   ((eta/2 + 1/2)*(zeta - 1)*(zeta + 1))/2,                                    -((eta/2 + 1/2)*(zeta - 1)*(zeta + 1))/2,                                       (zeta*(eta/2 - 1/2)*(zeta/2 + 1/2))/2,                                        -(zeta*(eta/2 - 1/2)*(zeta/2 + 1/2))/2,                                        -(zeta*(eta/2 + 1/2)*(zeta/2 + 1/2))/2,                                       (zeta*(eta/2 + 1/2)*(zeta/2 + 1/2))/2;...
                                       (zeta*(xi/2 - 1/2)*(zeta/2 - 1/2))/2,                                        -(zeta*(xi/2 + 1/2)*(zeta/2 - 1/2))/2,                                        -(zeta*(xi/2 - 1/2)*(zeta/2 - 1/2))/2,                                       (zeta*(xi/2 + 1/2)*(zeta/2 - 1/2))/2,                                    -((xi/2 - 1/2)*(zeta - 1)*(zeta + 1))/2,                                   ((xi/2 + 1/2)*(zeta - 1)*(zeta + 1))/2,                                   ((xi/2 - 1/2)*(zeta - 1)*(zeta + 1))/2,                                    -((xi/2 + 1/2)*(zeta - 1)*(zeta + 1))/2,                                       (zeta*(xi/2 - 1/2)*(zeta/2 + 1/2))/2,                                        -(zeta*(xi/2 + 1/2)*(zeta/2 + 1/2))/2,                                        -(zeta*(xi/2 - 1/2)*(zeta/2 + 1/2))/2,                                       (zeta*(xi/2 + 1/2)*(zeta/2 + 1/2))/2;...
 (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2) + (zeta*(xi/2 - 1/2)*(eta/2 - 1/2))/2, - (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 - 1/2) - (zeta*(xi/2 + 1/2)*(eta/2 - 1/2))/2, - (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2) - (zeta*(xi/2 - 1/2)*(eta/2 + 1/2))/2, (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 - 1/2) + (zeta*(xi/2 + 1/2)*(eta/2 + 1/2))/2, - (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta - 1) - (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta + 1), (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta - 1) + (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta + 1), (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta - 1) + (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta + 1), - (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta - 1) - (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta + 1), (xi/2 - 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2) + (zeta*(xi/2 - 1/2)*(eta/2 - 1/2))/2, - (xi/2 + 1/2)*(eta/2 - 1/2)*(zeta/2 + 1/2) - (zeta*(xi/2 + 1/2)*(eta/2 - 1/2))/2, - (xi/2 - 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2) - (zeta*(xi/2 - 1/2)*(eta/2 + 1/2))/2, (xi/2 + 1/2)*(eta/2 + 1/2)*(zeta/2 + 1/2) + (zeta*(xi/2 + 1/2)*(eta/2 + 1/2))/2];
else
    error('No shape functtion with thar order, sad panda');
end

end
