function D = hooke_trans_iso( EL, ET, nuLT, GLT, GTT )
%hooke_trans_iso Transverlsy isotropic elasticity
%   
nuTL = ET/EL*nuLT;

% Complience matrix
C = [    1/EL, -nuLT/EL, -nuLT/EL,  0    0    0
     -nuLT/EL,     1/ET, -nuTL/ET,  0    0    0
     -nuLT/EL, -nuTL/ET,     1/ET,  0    0    0
         0,        0,       0,    1/GLT, 0,   0
         0,        0,       0,      0, 1/GLT, 0
         0         0        0       0    0  1/GTT ];
 
D = inv(C); 

end

