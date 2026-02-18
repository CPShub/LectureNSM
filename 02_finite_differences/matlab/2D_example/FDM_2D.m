%% Code example 2: Finite Differences Method

% Lecture 'Numerical Simulation Methods'
% Felix Rutsch, May 2025

% solve stationary advection-diffusion equation with FDM

clearvars

%--- PARAMETERS (SI UNITS) ---
Lx = 1.0;              % [m]
Ly = 1.0;              % [m]
Nx = 11; Ny = 11;      % grid points
hx = Lx/(Nx-1);        % [m]
hy = Ly/(Ny-1);        % [m]

rho = 1.0e3;           % [kg/m^3], e.g. water
cp  = 4.18e3;          % [J/(kg K)], heat capacity
k   = 0.6;             % [W/(m K)], thermal conductivity
alpha = k/(rho*cp);    % [m^2/s], thermal diffusivity

v1 = 1.0; v2 = 0.5;   % [m/s]
% v1 = 0; v2 = 0;

%--- SPATIAL GRID ---
x = linspace(0, Lx, Nx);
y = linspace(0, Ly, Ny);
[X, Y] = meshgrid(x, y);

Q0 = 1000;
src_size = 1; % size parameter for the heat source

f = zeros(Nx,Ny);
cx = ceil(Nx/2);
cy = ceil(Ny/2);
f( cx-src_size:cx+src_size, cy-src_size:cy+src_size ) = Q0;

%--- FLATTEN SOURCE ---
Fv = reshape(f, Nx*Ny, 1);

%--- ASSEMBLE A and b ---
N = Nx*Ny;
A = zeros(N,N);
b = zeros(N,1);
idx = @(i,j) (j-1)*Nx + i;

for j = 2:Ny-1
  for i = 2:Nx-1
    p = idx(i,j);
    adv_x = rho*cp * v1 / (2*hx);
    adv_y = rho*cp * v2 / (2*hy);
    dfs_x = k      / hx^2;
    dfs_y = k      / hy^2;
    % central‐difference advection + diffusion
    A(p, idx(i,  j  )) = 2*(dfs_x+dfs_y);
    A(p, idx(i-1,j  )) = -dfs_x - adv_x;
    A(p, idx(i+1,j  )) = -dfs_x + adv_x;
    A(p, idx(i  ,j-1)) = -dfs_y - adv_y;
    A(p, idx(i  ,j+1)) = -dfs_y + adv_y;
    b(p) = Fv(p);
  end
end

%--- DIRICHLET BC: phi=300 K on all boundaries ---
T0 = 300;
for i = 1:Nx
  for j = [1 Ny]
    p = idx(i,j);
    A(p,p) = 1; b(p) = T0;
  end
end
for j = 2:Ny-1
  for i = [1 Nx]
    p = idx(i,j);
    A(p,p) = 1; b(p) = T0;
  end
end

%--- SOLVE & RESHAPE ---
phi_vec = A \ b;
phi = reshape(phi_vec, [Nx, Ny]); 

%--- VISUALIZE ---

% plot heat source
surf(X, Y, f');
xlabel('x'); ylabel('y'); zlabel('f(x,y)');
title('Heat source');

% surface plot
figure;
surf(X, Y, phi');
xlabel('x'); ylabel('y'); zlabel('\phi(x,y)');
title('2D Advection-Diffusion Solution');
shading interp;
colorbar;

% with grid
figure;
surf(X, Y, phi','EdgeColor','k');
view(2);
xlabel('x [m]'); ylabel('y [m]');
title('Temperature \phi(x,y) [K]');
shading interp;
colorbar;
axis equal tight;
