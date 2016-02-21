clear variables;
close all

problem = 'Konsol';
% problem = 'KonsolMedUtbredd';
%problem = 'InspandPlatta';

[mesh, elprop, M, bc, f ] = setup_problem(problem);
    
% Assemble
n = mesh.nel*(24)^2;
rows = zeros(n,1);
cols = zeros(n,1);
data = zeros(n,1);

nPassed = 1;
% TODO: Rewrite as a nonlinear-problem with Newton iterations /JB
for elIndex = 1:mesh.nel
    
    el(elIndex) = SolidShellLayered(3,3,5, mesh.ex(:,elIndex)', ...
        mesh.ey(:,elIndex)', mesh.ez(:,elIndex)', [2 2 3,2,3,3], M, elprop);
    [Ke, fe] = el(elIndex).computeLinearizedSystem([0 0 0]', [0 0 0]', elprop);
    % norm(Ke-Ke')/(norm(Ke)) -> Ke is symmetric
    
    % Assemble
    elDofs = mesh.edof(:,elIndex);
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


%Create stiffness matrix
K = sparse(rows,cols,data);

%Solve
a = solveq(K,f,bc);

%Element disp
ed = a(mesh.edof);
sfac = 1;
 exd = mesh.ex + ed(1:3:end,:)*sfac;
 eyd = mesh.ey + ed(2:3:end,:)*sfac;
 ezd = mesh.ez + ed(3:3:end,:)*sfac;
 
 figure(3);
 solid8draw(exd,eyd,ezd); hold on;
 view(3)
 axis equal
