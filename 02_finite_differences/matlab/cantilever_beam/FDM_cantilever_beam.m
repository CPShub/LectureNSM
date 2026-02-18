% Finite Difference Method

% Felix Rutsch, February 2025

% cantilever beam

clearvars;

%% parameters
N = 11;             % number of points
L = 1;              % [m] length of cantilever
q = 0;              % [N/m] line load
F = 1;              % [N] load at tip of cantilever
E = 1;              % [N/m2] Youngs modulus
I_const = 1;        % [m4] moment of inertia

delta_x = L/(N-1);
x_axis = 0:delta_x:L;

I = I_const*ones(N,1);   


%% solve

[Q, M, w] = FDM_solve_forward(N,L,q,F,E,I);


%% analytical solution
% for constant EI, F = 1, q(x) = 0
w_analytical = @(x) (F*L*x.^2/2 - F*x.^3/6)/(E*I_const);

%% plotting
plot_result(x_axis,Q,M,w,w_analytical)



