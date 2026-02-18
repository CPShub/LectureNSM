function [M, w] = FDM_solve_2nd_order(N,L,q,F,E,I)

delta_x = L/(N-1);

M_vector = -q*ones(N,1);
M_vector(N-1) = F; % BC: Q(x=L) = F
M_vector(N) = 0; % BC: M(x=L) = 0

M_coefs = sparse(N,N);
for i = 1:N-2
    M_coefs(i,i) = 1/(delta_x^2);
    M_coefs(i,i+1) = -2/(delta_x^2);
    M_coefs(i,i+2) = 1/(delta_x^2);
end
M_coefs(N,N) = 1; % BC M(x=L)
M_coefs(N-1,N-2:N) = [1, -4, 3]/(2*delta_x); % BC: Q(x=L) -> 2nd order backward difference

M = M_coefs\M_vector;


w_vector = -M./(E*I');
w_vector(1) = 0; % BC: w(x=0) = 0
w_vector(2) = 0; % BC: w'(x=0) = 0

w_coefs = sparse(N,N);
for i = 3:N
    w_coefs(i,i-2) = 1/(delta_x^2);
    w_coefs(i,i-1) = -2/(delta_x^2);
    w_coefs(i,i) = 1/(delta_x^2);
end
w_coefs(1,1) = 1; % BC w(x=0)
w_coefs(2,1:3) = [-3, 4, -1]/(2*delta_x); % BC: w'(x=0) -> 2nd order forward difference

w = w_coefs\w_vector;
