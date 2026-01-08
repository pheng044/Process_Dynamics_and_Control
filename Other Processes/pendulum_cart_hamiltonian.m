% ------------------------------------------------------------------------
% Patrick Heng
% 06/2025
% Script to simulate and plot an upright pendulum on a cart using
% Hamiltonian dynamics.
% ------------------------------------------------------------------------
close all; clear variables; clc


% Determine which plot to use, 1 = x,y plot, 2 = phase space
plot_flag = 1;

% ---- System Parameters -----
% Pendulum mass
m = 1;
% Cart mass
M = 1;
% Pendulum length
L = 1;
% Spring constant
k = 2;
% Gravitational constant
g = 9.81;

% ----- Simulation Parameters -----
% Minimum timestep
dt = 0.05;
% Initial condition x0 = [x,x_dot,theta,theta_dot]_0
x0 = [1,0,pi,0];



% ----- Solve dynamic model -----
f = @(t,x) pendulum_cart(x(1),x(2),x(3),x(4),m,M,L,k,g);

% Integrate ODE system
[t,x] = ode45(f,0:dt:100,x0);

% Pendulum x,y coordinates
xx = x(:,1) + L*sin(x(:,3));
yy = L*cos(x(:,3));

% Transform state variables to conjugate positions/momenta
H = hamiltonian(x(:,1),x(:,2),x(:,3),x(:,4),m,M,L,k,g);

% ----- Solution plotting -----
if plot_flag == 1
    % Animate pendulum in x,y plot
    % Physical limits to plot over
    x_lim = [min([xx,x(:,1)],[],'all'),max([xx,x(:,1)],[],'all')];
    y_lim = [min([yy;0],[],'all'),max([yy;0],[],'all')-0.9];

    % Loop through animation
    for i = 1:length(t)
        hold off
        % x,y trajectory
        plot(xx(1:i),yy(1:i),linewidth=1.5,color='b')
        hold on
        % Pendulum rod
        plot([x(i,1),xx(i)],[0,yy(i)],color='k',linewidth=2)
        % Pendulum mass
        plot(xx(i),yy(i),marker='.',color='r',markersize=30)
        % Cart center of mass
        plot(x(i,1),0*yy(i),marker='.',color='g',markersize=30)
        
        % Pretty plot parameters
        xlim(x_lim); ylim(y_lim)
        xlabel('$x$',interpreter='latex')
        ylabel('$y$',interpreter='latex')
        fontname('Serif'); fontsize(16,'points')
        grid on
        box on
        pause(0)
    end
else

    % Animate pendulum in phase space
    % Choose which cooridates to use for phase space
    % 1 = Cart position, 2 = cart momentum, 
    % 3 = pendulum position, 4 = pendulum momentum
    x1 = 2;
    x2 = 4;
    
    % Choose appropriate labels
    labels = {'$x$','$p_x$','$\theta$','$p_{\theta}$'};
    x_lab = labels(x1);
    y_lab = labels(x2);
    
    % Physical limits for plotting
    x_lim = [min(x(:,x1)),max(x(:,x1))];
    y_lim = [min(x(:,x2)),max(x(:,x2))];

    for i = 1:length(t)
        % Animate pendulum in phase space plot
        hold off
        % Trajectory up to time, t
        plot(x(1:i,x1),x(1:i,x2),linewidth=1.5,color='b')
        hold on
        % Trajectory terminating point
        plot(x(i,x1),x(i,x2),marker='.',color='r',markersize=30)

        % Pretty plot parameters
        xlim(x_lim); ylim(y_lim)
        xlabel(x_lab,interpreter='latex')
        ylabel(y_lab,interpreter='latex')
        fontname('Serif'); fontsize(16,'points')
        grid on
        box on
        pause(0)
    end
end


% Function to compute state space dynamics of pendulum on a cart
function dynamics = pendulum_cart(x,p_x,the,p_the,m,M,L,k,g)
    sin_the = sin(the);     % Precalculate sin(theta)
    cos_the = cos(the);     % Precalculate cos(theta)
    det = m*M*L^2 + m^2*L^2*sin_the.^2;     % Mass matrix determinant
    
    % Compute dynamics
    dynamics = [0;0;0;0];
    % Cart position
    dynamics(1) = (m*L^2*p_x - m*L*cos_the*p_the)*(det^-1);
    % Cart momentum
    dynamics(2) = -k*x;
    % Pendulum position
    dynamics(3) = (-m*L*cos_the*p_x + (m+M)*p_the)*(det^-1);
    % Pendulum momentum
    dynamics(4) = 0.5*m^2*L^2*sin_the*(m*L^2*cos_the*p_x^2 ...
                  - 2*(m+M)*L*p_x*p_the ...
                  + (m+M)*cos_the*p_the^2)*(det^-2) + m*g*L*sin_the;
    
end

% Function to compute Hamiltonian given conjugate coordinates
function H = hamiltonian(x,p_x,the,p_the,m,M,L,k,g)
    sin_the = sin(the);     % Precalculate sin(theta)
    cos_the = cos(the);     % Precalculate cos(theta)
    det = m*M*L^2 + m^2*L^2*sin_the.^2;     % Mass matrix determinant
    %mass = [m*L^2./det, -m*L*cos_the./det; -m*L*cos_the./det, (m+M)./det];

    % Generalized momentum vector
    p = [p_x;p_the];    
    % Compute kinetic energy
    T = 0.5*((m*L^2*det.^-1).*p_x.^2 ...
        - 2*(m*L*cos_the.*(det.^-1)).*p_x.*p_the ...
        + p_the.^2.*(m+M).*(det.^-1));
    % Return Hamiltonian
    H = m*g*L*cos_the + 0.5*k*x.^2 + T;

end