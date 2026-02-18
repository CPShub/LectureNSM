function [Q, M, w] = FDM_solve_forward(N,L,q,F,E,I)

delta_x = L/(N-1);

% coefficient matrix for BCs at i=N
coef_mtrx_1 = sparse(N,N);
for i=1:N-1
    coef_mtrx_1(i,i) = -1;
    coef_mtrx_1(i,i+1) = 1;
end
coef_mtrx_1(N,N) = 1; % BC: Q(x=L) and M(x=L)

%coefficient matrix for BCs at i=1
coef_mtrx_2 = sparse(N,N);
for i=2:N
    coef_mtrx_2(i,i-1) = -1;
    coef_mtrx_2(i,i) = 1;
end
coef_mtrx_2(1,1) = 1; % BC


% shear force
Q_vector = -q*delta_x*ones(N,1);
Q_vector(N) = F; % BC: Q(x=L) = F

Q = coef_mtrx_1\Q_vector;

% bending moment
M_vector = Q*delta_x;
M_vector(N) = 0; % BC: M(x=L) = 0

M = coef_mtrx_1\M_vector;

% incline
dw_vector = (-M*delta_x)./(E*I');
dw_vector(1) = 0;

dw = coef_mtrx_2\dw_vector;

% displacement
w_vector = dw*delta_x;
w_vector(1) = 0;

w = coef_mtrx_2\w_vector;
