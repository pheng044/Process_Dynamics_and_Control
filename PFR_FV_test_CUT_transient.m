% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 3 Mar 2025 (Steady State Solver)
% 14 Nov 2025 - 6 Dec 2025 (Transient Solver)
% Script to solve the PFR boundary value problem using a non-linear 
% finite volume scheme.
% Uses the CUT loop to solve the coupled transport problems:
%   C - Concentration equation
%   U - U momentum equation
%   T - Temperature equation
% which is iteratively solved until convergence
% -------------------------------------------------------------------

close all; clear all; clc;

% ----- INPUTS -----

% Initial concentration of A, cA0 (mol/m^3)
cA0 = 0.9*18.8;

% Initial concentrations of other species (mol/m^3)
% i = 1        2       3        4
%     Acetone  Ketene  Methane  Nitrogen (inert)
c0 = [cA0;0;0;0*cA0];

% Initial temperature, T0 (K)
T0 = 1000;

% Heat exchanger temperature, Ta (K)
Ta = 1150;

% Reactor length, L (m)
L = 4;

% Total Reactor Diameter (m)
D_tot = 1;

% Individual Tube Diameter (m)
D_tube = 0.0266446; % (1 inch, schedule 40 pipe)

% Number of tubes
N_tubes = 1000;

% Pipe roughness (m)
roughness = 0.00005;

% PFR velocity, U (m/s)
U0 = 8;

% Stoichiometric coefficients [a,b,c,d]
stoich = [-1;1;1;0];

% Molecular weights (kg/mol), 
% components arranged in same order as stoich vector
MW = [58.08;42.04;16.04;28.01]/1000;

% Heat capacity coefficients: Cp = a1 + a2*T + a3T^2 (J/mol*K), 
% components arranged in same order as stoich vector

a1 = [26.63;20.04;13.39;6.25];
a2 = [0.183;0.0945;0.077;8.78*1e-3];
a3 = [-45.86;-30.95;-18.71;-0.021]*1e-6;

% Arrhenius frequency factor, k0, (1/s)
k0 = exp(34.34);

% Activation energy, Ea, (J/mol)
Ea = 284500;

% Gas constant, R, (J/mol*K)
R = 8.314;

% Reference temperature, Tref (K)
Tref = 298;

% Enthalpy of reaction at Tref, Href (J/mol*K)
H_f_298 = [-216670;-61090;-71840;0];
Href = stoich'*H_f_298;

% Mixture average thermal conductivity, kappa (J/m*K)
kappa = 2e-1;

% Mixture average diffusivity, D (m^2/s)
D = 1e-2;

% Effective heat transfer coefficient, h = UA/V where U is the overall
% heat transfer coefficient (W/m^2*K), A is reactor surface area (m^2), 
% and V is the reactor volume (m^3)
h = 110*(4/D_tube);

% Mixture average viscocity, mu (kg/m*s)
mu = 30*1e-6;

% +++++ REVISION START (11/29/25) ++++++++++++++++++++++++++++++++++++
% Density of heating fluid (kg/m^3)
rho_a = 1.17469;
% Heat capacity of heating fluid (J/kg-K) 
Cp_a = 1159;

% Controller gain
Kc = 1000;
% Penalty coefficient to make T0 close to exit Ta
lambda = 1;
% Max mass flow rate for heating fluid (kg/s)
m_max = 1e1;

% Step input for CA
CA_input = @(t) -0.3*cA0*heaviside(t-1);
% +++++ REVISION END (11/29/25) +++++++++++++++++++++++++++++++++++++


% -------------------------------------------------------------------
% ----- SOLVER -----

% --- SOLVER PARAMETERS ---
% Number of nodes to solve for
nodes = 101;

% +++++ REVISION START (11/29/25) ++++++++++++++++++++++++++++++++++++
% Flag for real-time transient plotting
plot_transient = true;

% Steady momentum equation - flag for whether to solve the steady or 
% transient momementum equation
SS_U = true;
% +++++ REVISION END (11/29/25) ++++++++++++++++++++++++++++++++++++++

% Error tolerance for relative residuals, generally, anything below 
% 1e-6 is not recommended as it leads to oscillations in the iterations
% without further convergence
rel_tol = 1e-6;

% Max number of iterations for each time step
max_iter = 25;

% Max number of time steps
max_time_iter = 400;  
max_inner_iter = 50;    

% Below this tolerance, the solver will switch to a Newton algorithm,
% set to 0 if only fixed point iterations are desired. The Newton
% algorithm is efficient for problems below 10000 unknowns, but quickly
% becomes too slow for any bigger problems.
% Generally, do not set above 1e-3 since the Newton solver will converge 
% to unphysical solutions without a good initial guess
newton_tol = 0;

% If true, averaging will be performed on the outputs to make the 
% functions look visually smoother; introduces a small amount of error.
% However, by the mean value theorem, this error vanishes as dz -> 0.
smoothed_output = false;     

% Relaxation factors, 0.6 is good for many cases with a 1e-6 rel_tol,
% but they can be any number between 0 and 1
RFCA = 0.6;           % CA relaxation factor
RFT = RFCA;           % theta relaxation factor
RFU = RFCA;           % U relaxation factor

% -------------------------------------------------------------------
% -------------------------------------------------------------------
% ----- NO USER INPUTS BEYOND THIS POINT !!! ------------------------
% -------------------------------------------------------------------
% -------------------------------------------------------------------
% ----- PRE-SOLVE SETUP -----

% --- CALCULATED QUANTITIES ---
rho_0 = sum(c0.*MW);            % Initial mixture density
THETA = c0/cA0;                 % Normalized initial concentrations
% Initial mixture heat capacity
K = cA0*(THETA'*a1+THETA'*a2*T0+THETA'*a3*T0^2);           
p0 = sum(c0)*R*T0;              % Initial ideal gas pressure of mixture
p_max = 1e6*p0;
A_f = 0.25*pi*(D_tot^2-N_tubes*D_tube^2);  % Area for heating fluid


% Calulate dimensionless numbers
Pe_M = U0*L/D;                   % Mass Peclet number
Pe_T = K*U0*L/kappa;             % Thermal Peclet number
Da = k0*L/U0;                    % First order Damkohler number
Re = rho_0*U0*L/mu;              % Reynolds number
Eu = 1/(rho_0*U0^2);             % Euler number

% Calculate recurring constants
beta = -Ea/(R*T0);               % Dimensionless activation energy
TCp_ref = stoich'*a1*Tref+(1/2)*stoich'*a2*Tref^2+(1/3)*stoich'*a3*Tref^3;

Gamma = cA0*(Href - TCp_ref)/(K*T0);  % Dimensionless reference enthalpy
N = h*L/(K*U0);                       % Number of transfer units
epsilon = -sum(stoich)*cA0/(sum(c0)*stoich(1));  % Change in moles of rxn

% +++++ REVISION START (11/29/25) ++++++++++++++++++++++++++++++++++++
% Pipe friction factor
Re_D = rho_0*D_tube*U0/mu; % Reynolds number based on pipe diameter

if Re_D < 2300          % Laminar friction factor
    fric = 64/Re_D;     
else                    % Turbulent friction factor (Colebrook-White)
    fric = fsolve(@(x) x^-0.5 + ...
        2*log10(roughness/(3.7*D_tube)+2.51/(Re_D*sqrt(x))), ...
        64/Re_D)*L/D_tube;
end 
% +++++ REVISION END (11/29/25) +++++++++++++++++++++++++++++++++++++

% Generate a uniform computational mesh
z = linspace(0,1,nodes);

% Spatial finite difference
dz = z(2) - z(1);

% +++++ REVISION START (11/29/25) ++++++++++++++++++++++++++++++++++++
dt = dz;            % Time step chosen based on convective time scale
time = 0;           % Variable to store current time
% +++++ REVISION END (11/29/25) ++++++++++++++++++++++++++++++++++++++

% Store parameters in structure for easy usage in functions
params = struct('cA0',cA0,'T0',T0,'Ta',Ta,'stoich',stoich, ...
                'rho_0',rho_0,'THETA',THETA,'K',K,'p0',p0,'MW',MW, ...
                'Pe_M',Pe_M,'Pe_T',Pe_T,'Da',Da,'Re',Re,'Eu',Eu, ...
                'beta',beta,'Gamma',Gamma,'N',N,'epsilon',epsilon, ...
                'a1',a1,'a2',a2,'a3',a3,'dz',dz,'dt',dt,'nodes',nodes, ...
                'p_max',p_max,'fric',fric,'h',h,'U0',U0,'L',L, ...
                'D_tube',D_tube,'D_tot',D_tot,'N_tubes',N_tubes, ...
                'A_f',A_f, 'Kc',Kc,'rho_a',rho_a,'Cp_a',Cp_a, ...
                'm_max',m_max,'lambda',lambda,'mu',mu);


% --- SOLUTION INITIALIZATION ---
% Maximum number of time steps
inner_loop = 1:max_inner_iter;

% Initial guess for CA, assume an exponential curve
%CA = exp(-1000*z)';
CA = exp(-1000*z)';
%CA(1) = 1;

% Initial guess for theta, assume constant
theta = ones(nodes,1);

% Initial guess for U, assume constant
%U = 1e-2*ones(nodes,1);
U = 1 + 0*exp(-1000*z)';

% Initital guess for density field
rho = ones(nodes,1);

% Intitial guess for heating fluid field
theta_a = ones(nodes,1)*T0/Ta;

% Predict Jacobian sparsity pattern for non-linear solver
if newton_tol ~= 0 
    % Main diagonals
    J1 = spdiags(ones(2*nodes,3),[-1,0,1],nodes,nodes);
    
    % Extra off diagonals due to backwards and forward differencing 
    % at the boundaries
    J1(1,3) = 1;                % Left
    J1(nodes,nodes-2) = 1;      % Right
    
    jac_sparse = sparse([J1,J1,J1; J1,J1,J1; J1,J1,J1]);  
    clear J1

    % Generate options structure for fsolve with Jacobian sparsity pattern
    options = optimoptions('fsolve', Algorithm='trust-region', ...
        JacobPattern=jac_sparse,FunctionTolerance=rel_tol);
end 

% Count the number of linear and non-linear iterations used
Linear_Solves = 0;
Non_Linear_Solves = 0;

% Variable to store whether to use fixed point or Newton iteration.
% Start with fixed point iteration to establish convergence
newton = false;

% Preallocate residual histories
err_CA = zeros(max_inner_iter,1);
err_U = zeros(max_inner_iter,1);
err_T = zeros(max_inner_iter,1);

% Initial transport field profiles
X = 1 - U.*CA;                          % Extent of reaction
AA = THETA'*a1 + THETA'*a2*T0*theta + THETA'*a3*T0^2*theta.^2  ...
        - (stoich'*a1+(stoich'*a2*T0)*theta+(stoich'*a3*T0^2) ...
            *theta.^2).*(X/stoich(1));
H_t = (cA0/K)*AA.*theta./U;             % Enthalpy
U_mom_t = rho.*U;                       % Momentum

tt = (L/U0)*(0:dt:dt*max_time_iter);


% z-t plane variables
[Z,T] = meshgrid(z,tt);
CA_copy = zeros(nodes,max_time_iter);
T_copy = zeros(nodes,max_time_iter);
Ta_copy = zeros(nodes,max_time_iter);
flue_flow = zeros(max_time_iter,1);

% Initialize heating mass flow as 0
m = 0;


% -------------------------------------------------------------------
% ----- ITERATIVE SOLVER -----
tic
for t = 0:max_time_iter
% +++++ REVISION START (11/29/25) ++++++++++++++++++++++++++++++++++++

    if plot_transient == true
        % Plot variables over domain at current time step 
        hold off
        yyaxis left
        X = 1 - U.*CA;
    
        plot(z,CA,linewidth=2, color='b')            % Plot X
        ylabel('$C_A$',interpreter='latex')
        set(gca,'YColor','b')
        ylim([0,1])
        
        % Display physical time
        text(0.5,0.1,['$t=$ ' num2str(time,3) ' s'],interpreter='latex')
        
        yyaxis right
        plot(z,theta,linewidth=2, color='r')        % Plot theta
        ylabel('$\theta$, $U$, $\rho$, $\theta_a$',interpreter='latex')
        set(gca,'YColor','k')
        
        hold on
        plot(z,U,linewidth=2, color=[1,180,40]/256,linestyle='-')   % Plot U
        plot(z,rho,linewidth=2, color=[1,0.7,0],linestyle='-')   % Plot rho
        plot(z,theta_a,linewidth=2, color='m',linestyle='-')   % Plot theta_a
        ylim([0.75,1.25])
    
        % Pretty plot parameters
        xlabel('$z$',interpreter='latex'); 
        legend('$C_A$','$\theta$','$U$','$\rho$', '$\theta_a$', ...
            interpreter='latex', location='southwest')
        
        grid on; box on
        fontname('Serif'); fontsize(16,'points');
        
        pause(0)
    end

    CA_copy(:,t+1) = CA*cA0;
    T_copy(:,t+1) = theta*T0;
    Ta_copy(:,t+1) = theta_a*Ta;
    

    % Save state variables at previous time step for Euler scheme
    CA_t = CA;              % Concentration
    U_mom_t = rho.*U;       % U momentum
    % Enthalpy
    AA = THETA'*a1 + THETA'*a2*T0*theta + THETA'*a3*T0^2*theta.^2  ...
        - (stoich'*a1+(stoich'*a2*T0)*theta+(stoich'*a3*T0^2) ...
            *theta.^2).*(X/stoich(1));
    H_t = (cA0/K)*AA.*theta./U;
    % Heating fluid temperature
    theta_at = theta_a;
% +++++ REVISION END (11/29/25) ++++++++++++++++++++++++++++++++++++++   
for i = inner_loop

    if newton == false      % Linearized fixed point iteration
        
        % --- CA SOLVER ---
        % Save previous iteration to current
        CA_old = CA;

        % Gather CA stencil coefficients using values from the previous
        % iteration, perform a linear solve for the updated CA
        dist = CA_input(time);
        [A,f] = stencil_coefficients_CA(CA,U,theta,CA_t,dist,params);
        CA = A\f;
        %[l,u] = ilu(A);
        %CA = gmres(A,f,[],rel_tol,[],l,u,CA);
    
        % Update CA using underrelaxation
        CA = (1-RFCA)*CA_old + RFCA*CA; 

        % --- U SOLVER ---
        % Predict X, P, rho
        X = 1 - U.*CA;
        %X(U<delta) = delta;
        % Make sure X is bounded between 0 and 1, turn on if having
        % convergence issues in X
        % X(X<0) = 0; X(X>1) = 1;  

        P = pressure(X,U,theta,params);
        rho = density(X,theta,P,params);
        %rho(P==p_max) = 1;
        % Make sure density is positive, turn on if having convergence
        % issues
        %rho(rho<0) = 1e-14;
        
        % Save previous iteration to current
        U_old = U;
        % Gather U stencil coefficients using values from the previous
        % iteration, perform a linear solve for the updated U
        [A,f] = stencil_coefficients_U(U,rho,P,U_mom_t,SS_U,params);
        U = A\f;
        %[l,u] = ilu(A);
        %U = gmres(A,f,[],rel_tol,[],l,u,U);

        % Update U using underrelaxation
        U = (1-RFU)*U_old + RFU*U;
        
        % --- THETA SOLVER ---
        % Update X,P field using newly calculated U 
        X = 1 - U.*CA;
        %X(U<delta) = delta;
        %X(X<0) = 1e-14; X(X>1) = 1;
        P = pressure(X,U,theta,params);
        % Save previous iteration to current
        theta_old = theta;

        % Gather theta stencil coefficients using values from the previous
        % iteration, perform a linear solve for the updated theta
        [A,f] = stencil_coefficients_theta(CA,X,U,P,theta,theta_a,H_t, ...
                                                            params);
        theta = A\f;
        %[l,u] = ilu(A);
        %theta = gmres(A,f,[],rel_tol,[],l,u,theta);

        % Update theta using underrelaxation
        theta = (1-RFT)*theta_old + RFT*theta;

        % If theta drops below 0 anywhere, reset it to a small number. 
        % Required to prevent the divergence of the first few iterations 
        % with a poor initial guess.
        theta(theta<0) = 1e-14;
% +++++ REVISION START (11/29/25) ++++++++++++++++++++++++++++++++++++
        % --- THETA_A SOLVER ---
        % Save previous iteration to current
        theta_a_old = theta_a;

        % Gather theta stencil coefficients using values from the previous
        % iteration, perform a linear solve for the updated theta
        [A,f,m] = stencil_coefficients_Ta(theta,theta_a,z,theta_at,CA,params);
        theta_a = A\f;
        %[l,u] = ilu(A);
        %theta = gmres(A,f,[],rel_tol,[],l,u,theta);

        % Update theta using underrelaxation
        theta_a = (1-RFT)*theta_a_old + RFT*theta_a;
% +++++ REVISION END (11/29/25) ++++++++++++++++++++++++++++++++++++++
        % Increment the number of linear solves (4 per iteration)
        Linear_Solves = Linear_Solves + 4;      

    else    % Non-linear Newton solver
        
        % Update old variables with the previous iteration
        CA_old = CA;
        U_old = U;
        theta_old = theta;

        % Generate left hand side of F(x) = 0
        F = @(y) non_linear_PFR_discretization(y,params);

        % Use the previous iteration as the initial guess for the Newton
        % solver, place into a long column vector
        x0 = vertcat(CA,U,theta);
        
        % Non-linear Newton solver with MATLAB's fsolve
        x = fsolve(F,x0,options);

        CA = x(1:nodes);                % Extract CA from x vector
        U = x(nodes+1:2*nodes);         % Extract U from x vector
        theta = x(2*nodes+1:3*nodes);   % Extract theta from x vector

        % Increment the number of non-linear solves
        Non_Linear_Solves = Non_Linear_Solves + 1;

        % Stop the solver if the number of non-linear solves exceeds 3
        if Non_Linear_Solves > 3
            warning(['More than 3 non-linear solves, the solution ' ...
                'may not be converging...'])
            break
        end
        
    end
    
    % Calculate L2 relative residuals
    res_CA = norm(CA-CA_old,2)/norm(CA,2);
    res_U = norm(U-U_old,2)/norm(U,2);
    res_T = norm(theta-theta_old,2)/norm(theta,2);

    % Store residual history
    err_CA(i) = res_CA; err_U(i) = res_U; err_T(i) = res_T;

    % If the residuals are below the desired tolerance, switch to the
    % non-linear Newton solver
    if res_CA < newton_tol && res_T < newton_tol
        newton = true;
    end

    % Convergence criterion: all residuals less than the relative 
    % tolerance
    if res_CA < rel_tol && res_U < rel_tol && res_T < rel_tol
        break
    end

end

% +++++ REVISION START (11/29/25) ++++++++++++++++++++++++++++++++++++
%dt = min(dz./U);       % Turn on for adaptive time step
%params.dt = dt;     

% Convert non-dimensional time to physical time
time = time + dt*L/U0;    
flue_flow(t+1) = m;

end
% +++++ REVISION END (11/29/25) ++++++++++++++++++++++++++++++++++++++

% -------------------------------------------------------------------
% ----- POST PROCESSING -----

% Return the correct sized residual histories
err_CA = nonzeros(err_CA);
err_U = nonzeros(err_U);
err_T = nonzeros(err_T);

% Perform Gaussian averaging on the internal nodes of the outputs, 
% leaves the boundary nodes as is
if smoothed_output == true
    DD = (1/4)*spdiags([ones(nodes,1), 2*ones(nodes,1), ones(nodes,1)],...
                                [-1,0,1],nodes,nodes);
    DD(1,:) = 0; DD(1,1) = 1; DD(nodes,:) = 0; DD(nodes,nodes) = 1;
    X = DD*X;
    U = DD*U;
    theta = DD*theta;
end

% Clear excess variables
%clear R kappa D Ea Href Pe_M Pe_T Da Re K N Ta rho_0 mu TCp_ref
%clear T0 dCp Cp Gamma F dt_CA dt_T a1 a2 a3 c0 h MW 

% --- PLOT X, THETA, U, P, IN NON-DIMENSIONALIZED COORDINATES ---
figure
yyaxis left
X = 1 - U.*CA;
plot(z,X,linewidth=2, color='b')            % Plot X
ylabel('$X$',interpreter='latex')
set(gca,'YColor','b')
ylim([0,1])

yyaxis right

plot(z,theta,linewidth=2, color='r')        % Plot theta
ylabel('$\theta,\ U$',interpreter='latex')
set(gca,'YColor','k')

hold on
plot(z,U,linewidth=2, color=[1,180,40]/256,linestyle='-')   % Plot U

% Pretty plot  parameters
xlabel('$z$',interpreter='latex'); 
legend('$X$','$\theta$','$U$', interpreter='latex', location='best')

grid on; box on
fontname('Serif'); fontsize(16,'points');


% --- PLOT CONCENTRATION PROFILES ---
figure
X(U<1e-2) = 0;
c_A = cA0*(1-X)./U;
c_B = cA0*(THETA(2)-stoich(2)*X/stoich(1))./U;
c_C = cA0*(THETA(3)-stoich(3)*X/stoich(1))./U;

hold on
x = L*z';

plot(x,c_A,linewidth=2,color='b')
plot(x,c_B,linewidth=2,color=[1,180,40]/256)
plot(x,c_C,linewidth=2,color='r')

% Pretty plot  parameters
xlabel('$x$ (m)',interpreter='latex'); 
ylabel('$c_A$, $c_B$, $c_C$ (mol/m$^3$)',interpreter='latex');
legend('$c_A$', '$c_B$', '$c_C$',interpreter='latex', location='best')

grid on; box on
fontname('Serif'); fontsize(16,'points');

% --- PLOT MOLE FRACTION PROFILES ---
figure
X(U<1e-2) = 0;
c_T = c_A + c_B + c_C;
x_A = c_A./c_T;
x_B = c_B./c_T;
x_C = c_C./c_T;

hold on
x = L*z';

plot(x,x_A,linewidth=2,color='b')
plot(x,x_B,linewidth=2,color=[1,180,40]/256)
plot(x,x_C,linewidth=2,color='r')

% Pretty plot  parameters
xlabel('$x$ (m)',interpreter='latex'); 
ylabel('$x_A$, $x_B$, $x_C$',interpreter='latex');
legend('$x_A$', '$x_B$', '$x_C$',interpreter='latex', location='best')

grid on; box on
fontname('Serif'); fontsize(16,'points');

% --- PLOT REACTION RATE PROFILE ---
figure
rate = k0*cA0*exp(beta./theta).*CA;

plot(x,rate,linewidth=2,color='b')

% Pretty plot  parameters
xlabel('$x$ (m)',interpreter='latex'); 
ylabel('$-r_A$ (mol/m$^3$/s)',interpreter='latex');
grid on; box on
fontname('Serif'); fontsize(16,'points');
% +++++ REVISION START (11/29/25) ++++++++++++++++++++++++++++++++++++
% --- PLOT HEAT EXCHANGE PROFILES ---
figure
plot(x,theta_a*Ta,'-b',x,theta*T0,'-r',linewidth=2)

% Pretty plot  parameters
xlabel('$x$ (m)',interpreter='latex'); 
ylabel('$T$ (K)',interpreter='latex');
legend('HX Fluid','Reactor Fluid',interpreter='latex')
grid on; box on
fontname('Serif'); fontsize(16,'points');
% +++++ REVISION END (11/29/25) ++++++++++++++++++++++++++++++++++++++

% --- PLOT RESIDUALS ---
figure
semilogy(err_CA,linewidth=2,color='b')      % CA error

hold on
semilogy(err_T,linewidth=2,color='r')       % Theta error

semilogy(err_U,linewidth=2,color=[1,180,40]/256)    % U error

% Plot order of convergence
k = size(err_T,1);
%semilogy(1:k,2.^-(1:k), linewidth=2, linestyle='--', color='k')

% Pretty plot parameters
set(gca,'YColor','k')
xlabel('Iterations, $k$',interpreter='latex')
ylabel('$r_X,\ r_\theta,\ r_U$',interpreter='latex')
legend('$r_X$', '$r_\theta$', '$r_U$', ...
    interpreter='latex', location='southwest')

grid on; box on
fontname('Serif'); fontsize(16,'points');

save('PFR_dynamic_data.mat','Z','T','flue_flow','CA_copy', ...
    'T_copy','Ta_copy')


toc