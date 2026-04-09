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
%Xe = [0 1 1 0;
%      0 0 1 1];
Xe = [0 0 1 1;
      0 1 1 0]; % node coordinates
Fe = [1 1 1 1] * 12;

kappa = 6;
qp = [-1  1 -1  1;
      -1 -1  1  1] / sqrt(3); % gauss points

qw = [ 1  1  1  1];  % ?    

[Ke, be] = elem_heat_Q4(Xe,Fe,kappa,qp,qw)

function [Ke, be] = elem_heat_Q4(Xe,Fe,kappa,qp,qw)

% Initalization
Ne = 4;                 % Size of Ke, be
Ke = zeros(Ne,Ne);      % Element stiffness matrix
be = zeros(Ne,1);       % Element force vector
qn = size(qp,2);        % Number of quadrature points
    
% Quadrature loop
for k = 1:qn

    % Quadrature point
    Xi = qp(:,k);
    
    % Evaluation of parametric gradient shape functions Ni
    dNdXi = 0.25 * [ -1+Xi(2), -1+Xi(1);
                      1-Xi(2), -1-Xi(1);
                      1+Xi(2),  1+Xi(1);
                     -1-Xi(2),  1-Xi(1) ];

    % Jacobian (parametric gradient of coordinate transformation)
    J = Xe * dNdXi;
    Jinv = inv(J);
    detJ = abs(det(J));

    % Physical gradient of shape functions Ni 
    dNdX = dNdXi * Jinv;     

    % Integrand evaluation for element stiffness matrix
    Kek = dNdX * kappa * dNdX' * detJ;
    Ke = Ke + qw(k) * Kek;

    
    % element load vector

    % Evaluation of shape functions Ni
    Neval = 0.25 * [ (1-Xi(1))*(1-Xi(2));
                     (1+Xi(1))*(1-Xi(2));
                     (1+Xi(1))*(1+Xi(2));
                     (1-Xi(1))*(1+Xi(2)) ];

    % Integrand evaluation for element force vector
    bek = Neval * (Fe * Neval) * detJ;
    be = be + qw(k) * bek;
end

end