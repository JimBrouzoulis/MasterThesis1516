%Clear stuff
clc, clear variables;

%Mesh parameters
lx = 0.1; ly=0.01; lz = 0.002; 
nelx = 10; nely=1; nelz=1;

%Force
P = -40; 

%Mesh
[edof,coord,ex,ey,ez,dof,nel,ndofs,nno, side1nodes, side2nodes,...
 side3nodes, side4nodes, side5nodes] = cubeMesher(lx,ly,lz,nelx,nely,nelz);
neldofs = 24;

%Material
E = 100e9; nu = 0;
D = hooke(4,E,nu);

%Create reference to solid shell
element = SolidShell(3,3,3);

%Assemble
n = nel*(neldofs)^2;
rows = zeros(n,1);
cols = zeros(n,1);
data = zeros(n,1);

nPassed = 1;
f=zeros(ndofs,1);
for elIndex = 1:nel
    element.computeR(
    [Ke,fe] = solid8Eas2(ex(:,elIndex)', ey(:,elIndex)', ez(:,elIndex)',D, M, [0,0,0]', [3,3,3]);
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












