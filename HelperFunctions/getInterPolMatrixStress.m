function [P] = getInterPolMatrixStress(opt)

if(opt == 1)
    P = @(xi,eta,zeta) [eye(2), [0;0], [eta,0;0,xi], [0, zeta 0 eta*zeta, 0;0,0,zeta,0, xi*zeta], zeros(2,8);...
                        zeros(1,10), 1 0 0, xi*eta, 0 0 0 0;...
                        0 0 1 0 0 zeta 0 0 0 0 0 0 0 0 0 0 0 0;...
                        zeros(2,11), eye(2), [0;0], [eta 0 xi^2*zeta^2 0; 0, xi, 0, eta^2*zeta^2]];
elseif(opt == 2)

end

end