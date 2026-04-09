% -------------------------------------------------------------------------
% Introduction to the Finite Element Method
% Prof. Dr. Oliver Weeger, TU Darmstadt
% -------------------------------------------------------------------------
%
% Evaluation of element stiffness matrix & force vector of 
% bilinear quad element (Q4) for the heat equation
%

%Xe = [0 1 2 1;
%      0 0 1 1];
Xe = [0 1 1 0;
      0 0 1 1];
Fe = [1 1 1 1] * 12;
kappa = 6;
qp = [-1  1 -1  1;
      -1 -1  1  1] / sqrt(3);
qw = [ 1  1  1  1];      

[Ke, be] = elem_heat_Q4(Xe,Fe,kappa,qp,qw)

function [Ke, be] = elem_heat_Q4(Xe,Fe,kappa,qp,qw)

% Initalization
Ne = 4;                 % Size of Ke, be
Ke = zeros(Ne,Ne);      % Element stiffness matrix
be = zeros(Ne,1);       % Element force vector
qn = size(qp,2);        % Number of quadrature points



for i = 1:qn
    
    

end




end