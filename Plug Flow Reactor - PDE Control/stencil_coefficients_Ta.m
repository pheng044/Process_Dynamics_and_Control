% -------------------------------------------------------------------
% Patrick Heng
% 11/20/25 - 12/09/25
% 1D stencil coefficients for non-dimensionalized PFR problem
%  -> Heating fluid equation (theta_a)
% -------------------------------------------------------------------

function [A,f,m] = stencil_coefficients_Ta(theta,theta_a,z,theta_at,CA,params)

    % Load parameters needed to create coefficients
    dt = params.dt;
    dz = params.dz;
    nodes = params.nodes;
    U0 = params.U0;
    T0 = params.T0;
    h = params.h;
    L = params.L;
    Ta = params.Ta;
    Cp = params.Cp_a;
    rho = params.rho_a;
    A = params.A_f;
    Kc = params.Kc;
    lambda = params.lambda;
    m_max = params.m_max;
    
    % Calculate mass flow rate based on control law
    m = min(Kc*trapz(z,(theta*T0/Ta-theta_a).^2 + lambda*(theta-1).^2),m_max)

    % Calculate coefficients for linear system
    k1 = m/(rho*A*U0);
    k2 = h*L/(rho*Cp*U0);
    Pe_T = params.Pe_T;

    % Flux limiter
    r = zeros(nodes,1);     % Flux limiter variable
    r(2:nodes-1) = (theta_a(2:nodes-1)-theta_a(1:nodes-2)) ...
                    ./(theta_a(3:nodes)-theta_a(2:nodes-1));
    phi = max(zeros(nodes,1),min(ones(nodes,1),r));
    phi_1 = 1 - phi; 

    % Compute stencil coefficients
    W =  phi.*k1/(2*dz) - (1/Pe_T)*(1/dz^2);
    C =  1/dt + phi_1.*k1/dz + k2 + (1/Pe_T)*(2/dz^2);
    E = -phi.*k1/(2*dz) - phi_1.*k1/dz - (1/Pe_T)*(1/dz^2);
    
    % Return coefficients in convenient form for using spdiags
    diag = [W,C,E];
    
    % Generate A at internal nodes
    A = spdiags(diag,[-1,0,1],nodes,nodes+1);
    A = A(:,1:nodes);

    % Clear boundary rows
    %A(1,:) = 0; 
    A(nodes,:) = 0;
    
    % Left flux boundary
    A(nodes,nodes) = 1;
    % Right flux boundary
    %A(1,1) = 1/dz;
    %A(1,2) = -1/dz;
  
    % Forcing function with appropriate boundary conditions
    f = theta_at/dt + k2*theta*T0/Ta;     % Internal
    %f(1) = 0;                             % Left
    f(nodes) = 1;                          % Right

end
