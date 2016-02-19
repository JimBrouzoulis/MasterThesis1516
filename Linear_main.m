clear variables;
% close all;
%Main code for for assembling and setting bc:s and such.

%Meshing
lx = 0.1; ly=0.01; lz = 0.002;
nelx = 10; nely=1; nelz=1;

lx = 0.6; ly=0.3; lz = 0.01;
nelx = 60; nely=30; nelz=1;

%Force (Is [N] or [N/m^2] depending on your loading case
P = -40;
P = -0.2e6;

%Generate mesh
[edof,coord,ex,ey,ez,dof,nel,ndofs,nno,side1nodes,side2nodes,side3nodes, side4nodes,side5nodes] = cubeMesherHigherOrder(lx,ly,lz,nelx,nely,nelz,2,2,2,3);
neldofs = 3*(2*2*2);

%Material
E = 100e9; nu = 0;
D = hooke(4,E,nu);

%Eas interpolation
Mhat = getInterPolMatrix(1);

%Assembling
n = nel*(neldofs)^2;
rows = zeros(n,1);
cols = zeros(n,1);
data = zeros(n,1);

nPassed = 1;
f=zeros(ndofs,1);
for elIndex = 1:nel
    % Compute element stiffness
%     el(elIndex) = SolidShell(3,3,15,3,Mhat);
    el(elIndex) = SolidShell2(3,3,5, [2 2 3,2,3,3], Mhat);
    
    [Ke,fe] = el(elIndex).computeLinearizedSystem(ex(:,elIndex)',ey(:,elIndex)',ez(:,elIndex)', [0,0,0]', [0 0 0]',D);
    elDofs = edof(:,elIndex);
    
    % Assemble
    for j = 1:neldofs
        for k = 1:neldofs
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
% [f, bc] = cubeBC( 'KonsolMedUtbredd', f, P*ly*lx, dof, side1nodes, side2nodes, side3nodes, side4nodes, side5nodes);
% [f, bc] = cubeBC( 'Konsol'            , f, P      , dof, side1nodes, side2nodes, side3nodes, side4nodes, side5nodes);
[f, bc] = cubeBC( 'InspandPlatta', f, P*lx*ly, dof, side1nodes, side2nodes, side3nodes, side4nodes, side5nodes);
%Create K
K = sparse(rows,cols,data);

%Solve
a = solveq(K,f,bc);

%Element disp.
ed = a(edof);











