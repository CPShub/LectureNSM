% -------------------------------------------------------------------------
% fem2d: Educational MATLAB code for finite element discretizations in 2D
% Cyber-Physical Simulation Group, TU Darmstadt
% Creator: Prof. Dr. Oliver Weeger
% -------------------------------------------------------------------------
%
% Plotting of FEA results on regular, structured mesh
% for bilinear quad element (Q4) for linear elasticity
%

function plot_mesh_Q4(X,u,lx,ufac,colors)

outX = reshape(X(1,:),lx+1,[]);
outY = reshape(X(2,:),lx+1,[]);
outUx = reshape(u(1:2:end),lx+1,[]);
outUy = reshape(u(2:2:end),lx+1,[]);
outUn = sqrt(outUx .* outUx + outUy .* outUy);
outUn = outUn / max(max(outUn));
outXp = outX + ufac*outUx;
outYp = outY + ufac*outUy;
outZ1 = ones(size(outX,1),size(outX,2));
outZ0 = zeros(size(outX,1),size(outX,2));

figure; hold on;
if (nargin > 4)
    colors = reshape(colors,lx+1,[]);
    mesh(outX,outY,outZ0,colors, ...
        'FaceColor', 'texturemap', 'LineWidth', 1.0, 'FaceAlpha', 0.6, ...
        'EdgeColor', 'black', 'EdgeAlpha', 0.6);
else
    mesh(outX,outY,outZ0, ...
        'FaceColor', 'none', 'LineWidth', 1.0, ...
        'EdgeColor', 'black', 'EdgeAlpha', 0.6);
end
mesh(outXp,outYp,outZ1,outUn, ...
    'FaceColor', 'interp', 'FaceAlpha', 0.9, ...
    'EdgeColor', 'black',  'EdgeAlpha', 0.9);
axis equal;
grid on;
caxis('auto')

end