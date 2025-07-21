%% Code example 4: Finite Elements Method

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
Nx = 50; Ny = 50;    % grid points
hx = Lx/Nx;            % [m]
hy = Ly/Ny;            % [m]

rho = 1000;              % [kg/m^3], mass density, e.g. 1.0e3 for water
cp  = 1000;              % [J/(kg K)], heat capacity, e.g. 4.18e3
k   = 1;               % [W/(m K)], thermal conductivity, e.g. 0.6 
%alpha = k/(rho*cp);   % [m^2/s], thermal diffusivity

vx = 0; vy = 0;    % [m/s]

method = 'UDS'; % for the convective term: UDS or CDS

%--- SPATIAL GRID ---
x = hx/2:hx:Lx-hx/2;  % cell‐center x
y = hy/2:hy:Lx-hy/2;  % cell‐center y

Q0 = 1000;
src_size = 0.2; % [m] size of the heat source
ncell = src_size/(2*hx);
cx = ceil(Nx/2);
cy = ceil(Ny/2);

f = zeros(Nx,Ny);
f( cx-ncell+1:cx+ncell, cy-ncell+1:cy+ncell ) = Q0;

%--- FLATTEN SOURCE ---
Fv = reshape(f, Nx*Ny, 1);

%--- ASSEMBLE A and b ---
N = Nx*Ny;
A = sparse(N,N);
b = zeros(N,1);
idx = @(i,j) (j-1)*Nx + i;

switch method
    case 'CDS' % central differencing scheme
    % CDS for both the advective and diffusive terms   
    dfs_x = k/hx^2;
    dfs_y = k/hy^2;
    adv_x = rho*cp * vx / (2*hx);
    adv_y = rho*cp * vy / (2*hy);

    aW = dfs_x + adv_x;
    aE = dfs_x - adv_x;
    aS = dfs_y + adv_y;
    aN = dfs_y - adv_y;

    % for Dirichlet BCs
    aW2 = 2*dfs_x + adv_x;
    aE2 = 2*dfs_x - adv_x;
    aS2 = 2*dfs_y + adv_y;
    aN2 = 2*dfs_y - adv_y;

    case 'UDS' % upwind differencing scheme
    % CDS for the diffusive term
    % UDS for the advective term (maximum operator)
    dfs_x = k/hx^2;
    dfs_y = k/hy^2;
    adv_x = rho*cp * vx / hx;
    adv_y = rho*cp * vy / hy;

    aW = dfs_x + max(0,adv_x);
    aE = dfs_x + max(0,-adv_x);
    aS = dfs_y + max(0,adv_y);
    aN = dfs_y + max(0,-adv_y);

    % for Dirichlet BCs
    aW2 = 2*dfs_x + max(0,adv_x);
    aE2 = 2*dfs_x + max(0,-adv_x);
    aS2 = 2*dfs_y + max(0,adv_y);
    aN2 = 2*dfs_y + max(0,-adv_y);
   
end

for j = 1:Ny
  for i = 1:Nx
    p = idx(i,j);
    
    if i > 1 
        A(p, idx(i-1, j)) = -aW;
        A(p, p) = A(p, p) + aW;
    end
    if i < Nx
        A(p, idx(i+1, j)) = -aE;
        A(p, p) = A(p, p) + aE;
    end
    if j > 1 
        A(p, idx(i, j-1)) = -aS;
        A(p, p) = A(p, p) + aS;
    end
    if j < Ny
        A(p, idx(i, j+1)) = -aN;
        A(p, p) = A(p, p) + aN;
    end
    b(p) = Fv(p);

  end
end


%--- BOUNDARY CONDITIONS ---
% D for Diriclet, N for Neumann
% (use U or any other letter for undefined)

% west
BC_W = 'D';
phi_W = 300;
fW = 0;

% east
BC_E = 'D';
phi_E = 300;
fE = 0;

% south
BC_S = 'D';
phi_S = 300;
fS = 0; 

% north
BC_N = 'D';
phi_N = 300;
fN = 0; 


% west boundary (i=1)
switch BC_W
case 'D'
for j = 1:Ny
    p = idx(1, j);
    A(p, p) = A(p, p) + aW2;           
    b(p) = b(p) + aW2*phi_W;
end
case 'N'
for j = 1:Ny
    p = idx(1, j);           
    b(p) = b(p) + fW/hx;
end
end

% east boundary (i=Nx)
switch BC_E
case 'D'
for j = 1:Ny
    p = idx(Nx, j);
    A(p, p) = A(p, p) + aE2;           
    b(p) = b(p) + aE2*phi_E;
end
case 'N'
for j = 1:Ny
    p = idx(1, j);           
    b(p) = b(p) + fE/hx;
end
end

% south boundary (j=1)
switch BC_S
case 'D'
for i = 2:Nx-1
    p = idx(i, 1);
    A(p, p) = A(p, p) + aS2;          
    b(p) = b(p) + aS2*phi_S;
end
case 'N'
for j = 1:Ny
    p = idx(1, j);           
    b(p) = b(p) + fS/hx;
end
end

% north boundary (j=Ny)
switch BC_N
case 'D'
for i = 2:Nx-1
    p = idx(i, Ny);
    A(p, p) = A(p, p) + aN2;           
    b(p) = b(p) + aN2*phi_N;
end
case 'N'
for j = 1:Ny
    p = idx(1, j);           
    b(p) = b(p) + fN/hx;
end
end


%--- SOLVE & RESHAPE ---
phi_vec = A \ b;
phi = reshape(phi_vec, [Ny, Nx])'; 

%--- VISUALIZE ---
visualize.heat_source(x,y,hx,hy,f)
visualize.temperature_interp(x,y,hx,hy,phi)
visualize.temperature(x,y,phi)
