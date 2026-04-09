function stationary_diffusion_2d()
    % Domain and mesh
    Lx = 1; Ly = 1;       % Domain size
    nx = 10; ny = 10;     % Number of elements in x and y
    [nodes, elements] = createMesh(Lx, Ly, nx, ny);

    % Diffusion coefficient and source
    D = 1.0;
    f = @(x, y) 1.0;  % Constant source term

    % Assemble system
    nNodes = size(nodes, 1);
    K = sparse(nNodes, nNodes);
    F = zeros(nNodes, 1);

    for e = 1:size(elements, 1)
        [ke, fe] = elementStiffnessAndLoad(nodes(elements(e, :), :), D, f);
        K(elements(e,:), elements(e,:)) = K(elements(e,:), elements(e,:)) + ke;
        F(elements(e,:)) = F(elements(e,:)) + fe;
    end

    % Apply Dirichlet BC (u = 0 on boundary)
    boundaryNodes = find( ...
        abs(nodes(:,1)) < 1e-12 | abs(nodes(:,1)-Lx) < 1e-12 | ...
        abs(nodes(:,2)) < 1e-12 | abs(nodes(:,2)-Ly) < 1e-12);
    u = zeros(nNodes,1);
    freeNodes = setdiff(1:nNodes, boundaryNodes);

    % Solve system
    u(freeNodes) = K(freeNodes, freeNodes) \ F(freeNodes);

    % Plot solution
    trisurf(delaunay(nodes(:,1), nodes(:,2)), nodes(:,1), nodes(:,2), u, 'EdgeColor', 'none');
    colorbar;
    title('Solution of 2D Diffusion Equation');
    xlabel('x'); ylabel('y'); zlabel('u');
    view(3);
end

function [nodes, elements] = createMesh(Lx, Ly, nx, ny)
    dx = Lx / nx;
    dy = Ly / ny;
    [x, y] = meshgrid(0:dx:Lx, 0:dy:Ly);
    nodes = [x(:), y(:)];

    elements = zeros(nx*ny, 4);
    count = 1;
    for j = 1:ny
        for i = 1:nx
            n1 = (j-1)*(nx+1) + i;
            n2 = n1 + 1;
            n3 = n2 + nx + 1;
            n4 = n1 + nx + 1;
            elements(count, :) = [n1 n2 n3 n4];
            count = count + 1;
        end
    end
end

function [ke, fe] = elementStiffnessAndLoad(coords, D, f)
    % 2x2 Gauss quadrature
    gp = [-1/sqrt(3), 1/sqrt(3)];
    w = [1, 1];

    ke = zeros(4, 4);
    fe = zeros(4, 1);

    for i = 1:2
        for j = 1:2
            xi = gp(i); eta = gp(j);
            [N, dNdxi] = shapeFunctions(xi, eta);
            J = dNdxi' .* coords;
            detJ = det(J);
            dNdxy = J \ dNdxi;

            ke = ke + D * (dNdxy' * dNdxy) * detJ * w(i) * w(j);
            x_y = N * coords; % physical coords
            fe = fe + N' * f(x_y(1), x_y(2)) * detJ * w(i) * w(j);
        end
    end
end

function [N, dNdxi] = shapeFunctions(xi, eta)
    % Shape functions for Q4 element
    N = 0.25 * [(1 - xi)*(1 - eta);
                (1 + xi)*(1 - eta);
                (1 + xi)*(1 + eta);
                (1 - xi)*(1 + eta)];

    dNdxi = 0.25 * [-(1 - eta), -(1 - xi);
                     (1 - eta), -(1 + xi);
                     (1 + eta),  (1 + xi);
                    -(1 + eta),  (1 - xi)]';
end