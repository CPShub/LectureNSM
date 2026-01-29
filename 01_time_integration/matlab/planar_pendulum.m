%% Code example 1: planar mathematical pendulum

% Lecture 'Numerical Simulation Methods'
% Felix Rutsch, May 2025

clearvars
close all

% Parameters
g = 9.81;                   % gravity [m/s^2]
L = 1.0;                    % length of pendulum [m]
delta_t = 0.05;             % time step [s]
T = 10;                     % total time [s]
N = round(T/delta_t);       % number of time steps

% Preallocate
phi = zeros(1, N+1);
v = zeros(1, N+1);
t = 0:delta_t:T;

% Initial conditions 
phi(1) = pi/4;              % [rad] e.g. pi/4 for 45 degrees
v(1) = 0;                   % [rad/s]

method = 'trapezoid rule';  % define the method
linearized = false;          % define usage of lin. or nonlin. pendulum eq.

switch method
case 'explicit Euler' %----------------------------------------------------
    for n=1:N
        phi(n+1) = phi(n) + v(n)*delta_t;

        if linearized
            v(n+1) = v(n) - g/L*phi(n)*delta_t;
        else
            v(n+1) = v(n) - g/L*sin(phi(n))*delta_t; 
        end
    end

case 'implicit Euler' %----------------------------------------------------
    for n = 1:N
        
        if linearized % linear system of equations        
            A = [1,-delta_t;
                 g/L*delta_t, 1];
            p = A \ [phi(n); v(n)];
            phi(n+1) = p(1);
            v(n+1) = p(2);
    
        else % non-linear system of equations 

            % Initial guess for Newton iteration
            phi_guess = phi(n);
            v_guess = v(n);
    
            % Newton iteration
            for iter = 1:10
                % Residuals
                R1 = phi_guess - phi(n) - delta_t * v_guess;
                R2 = v_guess - v(n) + delta_t * (g/L) * sin(phi_guess);
    
                % Jacobian
                J = [1,        -delta_t;
                     delta_t*(g/L)*cos(phi_guess), 1];
    
                % Update
                h = -J \ [R1; R2];
                phi_guess = phi_guess + h(1);
                v_guess = v_guess + h(2);
    
                if norm(h) < 1e-10
                    break;
                end
            end
    
            % Store updated values
            phi(n+1) = phi_guess;
            v(n+1) = v_guess;

        end
    end 

case 'symplectic Euler' %--------------------------------------------------
for n=1:N

    phi(n+1) = phi(n) + v(n)*delta_t;

    if linearized
        v(n+1) = v(n) - g/L*phi(n+1)*delta_t;
    else
        v(n+1) = v(n) - g/L*sin(phi(n+1))*delta_t; 
    end
end

case 'trapezoid rule' %----------------------------------------------------
for n = 1:N

    if linearized % linear system of equations        
            A = [1, -0.5*delta_t;
                 0.5*delta_t*(g/L), 1];
            B = [phi(n) + delta_t*v(n)/2; v(n) - delta_t*g/(2*L)*phi(n)];
            p = A \ B;
            phi(n+1) = p(1);
            v(n+1) = p(2);
    
    else % non-linear system of equations
    
        % Initial guess for Newton iterations
        phi_guess = phi(n);
        v_guess = v(n);
        
        for iter = 1:10
            R1 = phi_guess - phi(n) - 0.5 * delta_t * (v(n) + v_guess);
            R2 = v_guess - v(n) + 0.5 * delta_t * (g/L) * (sin(phi(n)) + sin(phi_guess));
    
            J = [1, -0.5*delta_t;
                 0.5*delta_t*(g/L)*cos(phi_guess), 1];
    
            h = -J \ [R1; R2];
            phi_guess = phi_guess + h(1);
            v_guess = v_guess + h(2);
    
            if norm(h) < 1e-10 
                break; 
            end
        end
        phi(n+1) = phi_guess;
        v(n+1) = v_guess;
    end
end


% -------------------------------------------------------------------------
otherwise
    error('unknown method name')

end

% Plot
visualize.angle(phi, v, t)

% energy
visualize.energy(phi, v, t, L, g) % only for nonlinear eq.

% compare with analytical solution
%visualize.analytical(phi,t,g,L) % only for linearized eq. and certain BCs

% animate
%visualize.animate(phi, t, L, delta_t)


%#ok<*UNRCH> % suppress warning unreachable code lines
