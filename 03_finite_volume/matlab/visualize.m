classdef visualize
   methods (Static)
        function temperature(x,y,phi)
            figure;
            surf(x,y,phi);

            shading interp;  % Enables smooth color transitions
            view(2);         % 2D top-down view

            axis equal tight;
            colormap(parula);  % color scheme, other options e.g. "hot", "jet"
            colorbar;

            xlabel('x'); ylabel('y'); zlabel('\phi(x,y)');
            title('Temperature Distribution');
        end

        function temperature_interp(x,y,hx,hy,phi)
            figure;
            imagesc(x,y,phi);
            set(gca,'YDir','normal') % correct y-axis direction
        
            hold on;
            % horizontal lines
            for row = 0:hx:1
            plot([0, 1], [row, row], 'k');
            end
            % vertical lines
            for col = 0:hy:1
                plot([col, col], [0, 1], 'k');
            end
            % optional: remove default axis ticks and add custom ones centered in cells
            % set(gca, 'XTick', 1:nCols, 'YTick', 1:nRows);
            hold off

            axis equal tight;
            colormap(parula);  % color scheme, other options e.g. "hot", "jet"
            colorbar;

            xlabel('x'); ylabel('y'); zlabel('\phi(x,y)');
            title('Temperature Distribution');
        end

        function heat_source(x,y,hx,hy,Q)
            figure;
            imagesc(x,y,Q);
            set(gca,'YDir','normal') % correct y-axis direction
        
            hold on;
            % horizontal lines
            for row = 0:hx:1
            plot([0, 1], [row, row], 'k');
            end
            % vertical lines
            for col = 0:hy:1
                plot([col, col], [0, 1], 'k');
            end
            % optional: remove default axis ticks and add custom ones centered in cells
            % set(gca, 'XTick', 1:nCols, 'YTick', 1:nRows);
            hold off

            axis equal tight;
            colormap(parula);  % color scheme, other options e.g. "hot", "jet"
            colorbar;

            xlabel('x'); ylabel('y'); zlabel('Q(x,y)');
            title('Heat Source');
        end
   end
end

