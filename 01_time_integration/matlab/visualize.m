classdef visualize
    
    methods (Static)
        function angle(phi, v, t)
            figure;
            subplot(2,1,1);
            plot(t, phi, 'b');
            ylabel('\phi (rad)');
            title('2D Mathematical Pendulum');

            subplot(2,1,2);
            plot(t, v, 'r');
            ylabel('v (rad/s)');
            xlabel('Time (s)');
        end

        function animate(phi, t, L, delta_t)
            figure;
            axis equal;
            axis([-1.2 1.2 -1.2 0.2] * L);
            hold on;
            grid on;
        
            % Elements to update
            blob = plot(0, 0, 'ro', 'MarkerSize', 12, 'MarkerFaceColor', 'r');
            rod = plot([0, 0], [0, 0], 'k-', 'LineWidth', 2);
            time_text = text(-1.1*L, 0.15*L, '', 'FontSize', 12);
        
            for n = 1:size(phi,2)
                x = L * sin(phi(n));
                y = -L * cos(phi(n));
        
                set(rod, 'XData', [0 x], 'YData', [0 y]);
                set(blob, 'XData', x, 'YData', y);
                set(time_text, 'String', sprintf('Time: %.2f s', t(n)));
        
                drawnow;
                pause(delta_t);
            end
        end

        function energy(phi, v, t, L, g)
            m = 1;  % assume unit mass
            kinetic = 0.5 * m * (L^2) * v.^2;
            potential = m * g * L * (1-cos(phi));
            total = kinetic + potential;

            figure;
            hold on
            plot(t, kinetic, 'LineWidth', 1.5);
            plot(t, potential, 'LineWidth', 1.5);
            plot(t, total, 'LineWidth', 1.5);
            xlabel('Time (s)');
            ylabel('Energy (J)');
            title('Energy');
            legend('kinetic energy','potential energy','total energy')
            grid on;
        end

        function analytical(phi,t,g,L)
            % implemented only for v(1) = 0, linearized equation
            omega = sqrt(g/L);
            analyt_solution = phi(1)*cos(omega*t); % initial condition defines the amplitude

            figure;
            subplot(2,1,1);
            hold on
            plot(t,analyt_solution,'b')
            plot(t, phi, 'r');
            ylabel('\phi (rad)');
            title('2D Mathematical Pendulum');
            legend('analytical solution','numerical solution')
            grid on

            subplot(2,1,2);
            title('2D Mathematical Pendulum');
            plot(t, abs(phi-analyt_solution), 'r');
            ylabel('numerical error (rad)');
            xlabel('Time (s)');

        end   
        
    end
end

