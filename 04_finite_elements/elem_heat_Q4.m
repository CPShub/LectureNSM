function [Ke, be] = elem_heat_Q4(Xe,Fe,kappa)

qp = [-1  1 -1  1;
      -1 -1  1  1] / sqrt(3); % gauss points

qw = [ 1  1  1  1]; % weigths of the gauss points

% Initalization
Ne = 4;                 % Size of Ke, be
Ke = zeros(Ne,Ne);      % Element stiffness matrix
be = zeros(Ne,1);       % Element force vector
qn = size(qp,2);        % Number of quadrature points
    
% Quadrature loop
for k = 1:qn

    % Quadrature point
    xi = qp(1,k);
    eta = qp(2,k);
    
    % Evaluation of parametric gradient shape functions Ni
    dNdXi = 0.25 * [ -1+eta, -1+xi;
                      1-eta, -1-xi;
                      1+eta,  1+xi;
                     -1-eta,  1-xi ];

    % Jacobian (parametric gradient of coordinate transformation)
    J = Xe * dNdXi;
    Jinv = inv(J);
    detJ = abs(det(J));

    % Physical gradient of shape functions Ni 
    dNdX = dNdXi * Jinv;     

    % Integrand evaluation for element stiffness matrix
    Kek = dNdX * kappa * dNdX' * detJ; % corresponds to weak form
    Ke = Ke + qw(k) * Kek; % multiply with weight of gauss point

    
    % element load vector

    % Evaluation of shape functions Ni
    Neval = 0.25 * [ (1-xi)*(1-eta);
                     (1+xi)*(1-eta);
                     (1+xi)*(1+eta);
                     (1-xi)*(1+eta) ];

    % Integrand evaluation for element force vector
    bek = Neval * (Fe * Neval) * detJ;
    be = be + qw(k) * bek;
end

end