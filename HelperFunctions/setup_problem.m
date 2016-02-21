function [m, elprop, M, bc, f ] = setup_problem(problem)
%setup problem: creates all the necessary data for a given problem 

% Corresponds to an 8 noded brick element with 3 dofs per node
nnoxel = 2; nnoyel = 2; nnozel = 2; ndofsno = 3; 


% Setup data for the problems
switch problem

    case 'InspandPlatta'
        
        % Mesh parameters
        lx = 0.1; ly=0.1; lz = 0.002;
        nelx = 10; nely=10; nelz=1;

        % Loading
        P = -0.2e6;
        load = P*lx*ly;
        
        %Material properties
        EL = 60E9;
        ET = 10E9;
        nuLT = 0.26;
        GLT = 5E9;
        GTT = 3.5E9;
        nuTL = ET/EL*nuLT;
        
        D = hooke_trans_iso( EL, ET, nuLT, GLT, GTT );
        D = hooke(4, 100e9, 0);
        angles = [0 0];
        
    case 'KonsolMedUtbredd'
        lx = 0.1; ly=0.01; lz = 0.002;
        nelx = 10; nely=1; nelz=1;
        load = 0;
        angles = [0];
        D = hooke(4,100e9,0);
        
    case 'Konsol'
        lx = 0.1; ly=0.01; lz = 0.002;
        nelx = 10; nely=1; nelz=1;
        load = -40;
        angles = [0];
        D = hooke(4,100e9,0);
        
    otherwise
        error('Unknown problem type')

end

% Create element properties
elprop = ElementProperties();
elprop.setup(lz, angles, D);

% Create mesh
m = Mesh();
m.create_cube_mesh(lx,ly,lz,nelx,nely,nelz,nnoxel,nnoyel,nnozel,ndofsno);

% Interpolation matrix
M = getInterPolMatrix(1);

% Load and boundary conditions
f=zeros(m.ndofs,1);

[f, bc] = cubeBC( problem, f, load, m.dof, m.side1nodes,...
    m.side2nodes, m.side3nodes, m.side4nodes, m.side5nodes);

end

