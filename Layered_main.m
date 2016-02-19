clear variables;
close all

%Meshing
lx = 0.6; ly=0.3; lz = 0.01;
nelx = 60; nely=30; nelz=1;

% lx = 0.1; ly=0.01; lz = 0.002;
% nelx = 10; nely=1; nelz=1;


%Force
P = -0.2e6;
% P = -40;

%Generate Mesh
[edof,coord,ex,ey,ez,dof,nel,ndofs,nno,side1nodes,side2nodes,side3nodes, side4nodes,side5nodes] = cubeMesherHigherOrder(lx,ly,lz,nelx,nely,nelz,2,2,2,3);
neldofs = 3*(2*2*2);

%Laminate properties
totalThickness = lz;        
angle = pi/180*[0 0]; 
nla = length(angle);
hlam = totalThickness/nla; % Thickness of a lamina
lamcoord =  0:hlam:totalThickness; % Coordinates for interfaces between lamina

%Properties for elementroutine
para.angles = angle;
para.coords = lamcoord;

%Material properties
EL = 60E9;    
ET = 10E9;   
nuLT =0.26;    
GLT = 5E9;  
GTT = 3.5E9;     
nuTL = ET/EL*nuLT;

%Interpolation matrix
M = getInterPolMatrix(1);

%Complience matrix
ComplianceMatrix =  [1/EL, -nuLT/EL, -nuLT/EL, 0 0 0;...
                     -nuLT/EL, 1/ET, -nuTL/ET, 0 0 0;...
                     -nuLT/EL, -nuTL/ET, 1/ET, 0 0 0;...
                     0,0,0, 1/GLT, 0, 0;...
                     0,0,0,0, 1/GLT,0;...
                     0 0 0 0 0 1/GTT];
 
D = inv(ComplianceMatrix); 
D = hooke(4,100e9,0);

para.D_LT = D;

%Assemble
n = nel*(24)^2;
rows = zeros(n,1);
cols = zeros(n,1);
data = zeros(n,1);

nPassed = 1;
f=zeros(ndofs,1);
for elIndex = 1:nel
    
    el(elIndex) = SolidShellLayered(3,3,5, ex(:,elIndex)', ey(:,elIndex)', ez(:,elIndex)', [2 2 3,2,3,3], M, para);
    [Ke, fe] = el(elIndex).computeLinearizedSystem([0 0 0]', [0 0 0]');
    
    % Assemble
    elDofs = edof(:,elIndex);
    for j = 1:24
        for k = 1:24
            rows(nPassed) = elDofs(j);
            cols(nPassed) = elDofs(k);
            data(nPassed) = Ke(j,k);
            nPassed = nPassed + 1;
        end
    end
    f(elDofs) = f(elDofs) + fe;
end

%Boundary condition
% [f, bc] = cubeBC( 'KonsolMedUtbredd', f, 0*P*ly*lx, dof, side1nodes, side2nodes, side3nodes, side4nodes, side5nodes);
[f, bc] = cubeBC( 'InspandPlatta', f, P*lx*ly, dof, side1nodes, side2nodes, side3nodes, side4nodes, side5nodes);
% [f, bc] = cubeBC( 'Konsol', f, P, dof, side1nodes, side2nodes, side3nodes, side4nodes, side5nodes);

%Create stiffness matrix
K = sparse(rows,cols,data);

%Solve
a = solveq(K,f,bc);

%Element disp
ed = a(edof);

% exd = ex + ed(1:3:end,:)*50;
% eyd = ey + ed(2:3:end,:)*50;
% ezd = ez + ed(3:3:end,:)*50;
% 
% figure(3);
% solid8draw(exd,eyd,ezd); hold on;
% view(3)
% axis equal
