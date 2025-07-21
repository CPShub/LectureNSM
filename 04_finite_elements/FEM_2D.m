%% Code example 4: Finite Element Method

% Lecture 'Numerical Simulation Methods'
% Felix Rutsch, July 2025

% solve stationary diffusion equation with FEM

% y
% ^
% |
% o ---> x

clearvars
close all

%--- PARAMETERS (SI UNITS) ---
Lx = 1.0;              % [m]
Ly = 1.0;              % [m]
nx = 10;               % elements in x direction
ny = 10;               % elements in y direction

k   = 1;               % [W/(m K)], thermal conductivity 
%alpha = k/(rho*cp);   % [m^2/s], thermal diffusivity

T0 = 300;              % [K], Temperature at the Dirichlet boundaries

Q0 = 1000;

% --- GENERATE FE MESH ---
[nodes, edof] = mesh_rect_Q4(Lx, Ly, nx, ny); % edof = "element degrees of freedom"

% elements for heat source
a = ny*nx/2 - ny/2;
Q0_elems = [a, a+1, a+ny, a+ny+1];

%--- ASSEMBLY ---
nnodes = (nx+1)*(ny+1);         % number of nodes
nelem = nx*ny;                  % number of elements

K = sparse(nnodes, nnodes);     % initialize stiffness matrix
F = zeros(nnodes, 1);           % initialize load vector

% assembly loop
for e = 1:nelem

    Xe = nodes(:,edof(e,:));

    if ismember(e, Q0_elems)
        Fe = [1 1 1 1]*Q0*(Lx/nx)*(Ly/ny)/4;
    else
        Fe = [0 0 0 0];
    end

    % element stiffness matrix
    [ke, fe] = elem_heat_Q4(Xe,Fe,k);

    % assembly
    K(edof(e,:), edof(e,:)) = K(edof(e,:), edof(e,:)) + ke;
    F(edof(e,:)) = F(edof(e,:)) + fe;
end

%--- BOUNDARY CONDITIONS ---
boundaryNodes = [];
for n=1:nnodes
    if nodes(1,n) == -Lx/2 || nodes(1,n) == Lx/2 || ...
       nodes(2,n) == -Ly/2 || nodes(2,n) == Ly/2          
    boundaryNodes = [boundaryNodes n];
    end
end
freeNodes = setdiff(1:nnodes, boundaryNodes);

phi = zeros(nnodes,1);
phi(boundaryNodes) = T0; % Dirichlet BC

% Modify load vector to account for Dirichlet BC
F = F - K * phi;

%--- SOLVE ---
phi(freeNodes) = K(freeNodes, freeNodes) \ F(freeNodes);

%--- VISUALIZE ---
visualize_FEM(nodes', edof, phi)

%u_matrix = reshape(u,[ny+1,nx+1]);

