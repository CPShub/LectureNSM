% convergence study

clearvars
close all

% Parameters
g = 9.81;                   % gravity [m/s^2]
L = 1.0;                    % length of pendulum [m]
%delta_t = 0.05;             % time step [s]
%T = 1;                     % total time [s]
%N = round(T/delta_t);       % number of time steps

stepsize = [0.5 0.1 0.05 0.01 0.005 0.001 0.0005 0.0001 0.00005 0.0001];

N = 2;
% Preallocate
phi = zeros(1, N+1);
v = zeros(1, N+1);

% Initial conditions 
phi(1) = pi/10;              % [rad] e.g. pi/4 for 45 degrees
v(1) = 0;                   % [rad/s]

for i=1:10

    delta_t = stepsize(i);
    t = 0:delta_t:2*delta_t;
    N = 2;

    % analytical solution
    % implemented only for v(1) = 0, linearized equation
    omega = sqrt(g/L);
    analyt_phi = phi(1)*cos(omega*t);% initial condition defines the amplitude
    analyt_v = -phi(1)*omega*sin(omega*t);

    % implicit Euler
    for n=1:N
        A = [1,-delta_t;
                 g/L*delta_t, 1];
            p = A \ [phi(n); v(n)];
            phi(n+1) = p(1);
            v(n+1) = p(2);
    end
    error_impl(i) = phi(2)-analyt_phi(2);


    % symplectic Euler
    for n=1:N
        phi(n+1) = phi(n) + v(n)*delta_t;
        v(n+1) = v(n) - g/L*phi(n+1)*delta_t;
    end
    error_sympl(i) = phi(2)-analyt_phi(2);
    
    % trapezoid rule
    for n = 1:N
        A = [1, -0.5*delta_t;
             0.5*delta_t*(g/L), 1];
        B = [phi(n) + delta_t*v(n)/2; v(n) - delta_t*g/(2*L)*phi(n)];
        p = A \ B;
        phi(n+1) = p(1);
        v(n+1) = p(2);
    end

    error_trpz(i) = phi(2)-analyt_phi(2);

end

loglog(stepsize,error_sympl,'-o','DisplayName','Symplektisches Euler-Verfahren')
hold on
loglog(stepsize,error_trpz,'-o','DisplayName','Trapezregel')
grid on
ylabel('Fehler |y*-y| im ersten Zeitschritt')
xlabel('Zeitschrittweite \Delta t')
legend('Location','southeast')

p1 = ( log(error_trpz(2)/error_trpz(4)) )/( log(stepsize(2)/stepsize(4)) )