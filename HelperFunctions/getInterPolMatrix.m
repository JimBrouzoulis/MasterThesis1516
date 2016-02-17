function [M] = getInterPolMatrix(opt)

if(opt == 1)
    M = @(xi,eta,zeta) [zeros(2,3); [0 0 zeta] ;zeros(1,3); [xi,0,0; 0,eta,0]];
% M = @(xi,eta,zeta) [zeros(4,1); [xi;0]];
elseif(opt == 2)
    M = @(xi,eta,zeta) [zeros(4,8); [xi,zeta,0 0 xi*eta, xi*zeta,0,0; 0 0 eta, zeta, 0 0 xi*eta, eta*zeta]];
elseif(opt == 3) 
    M = @(xi,eta,zeta) [zeros(4,12); [xi 0 xi*eta 0 (1/5 - zeta^2) 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2), 0 0 xi*eta*(1/5-zeta^2), 0; 0 eta 0 xi*eta 0 (1/5 - zeta^2) 0 0 xi*(1/5-zeta^2), eta*(1/5-zeta^2) 0 xi*eta*(1/5-zeta^2)]];  
%     M = @(xi,eta,zeta) [zeros(4,12); [xi 0 xi*eta 0 (-2/5 - zeta^2) 0 xi*(-2/5-zeta^2), eta*(-2/5-zeta^2), 0 0 xi*eta*(1/5-zeta^2), 0; 0 eta 0 xi*eta 0 (-2/5 - zeta^2) 0 0 xi*(-2/5-zeta^2), eta*(-2/5-zeta^2) 0 xi*eta*(-2/5-zeta^2)]];  
elseif(opt == 4)
    M = @(xi,eta,zeta) [zeros(2,4); zeta, 0, 0, 0; zeros(1,4); [xi,0,xi*eta,0; 0,eta,0,xi*eta]];
elseif(opt == 5)
     M = @(xi,eta,zeta) [zeros(4,6); [1 0, xi,0,xi*eta,0; 0, 1, 0,eta,0,xi*eta]];
% M = @(xi,eta,zeta) [0,0,0,0,xi,0]';
end

end