% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025 (Steady State Solver)
% 14 Nov 2025 - 6 Dec 2025 (Transient Solver)
% 1D stencil coefficients for non-dimensionalized PFR problem
%  -> theta equation
% -------------------------------------------------------------------

function [A,f] = stencil_coefficients_theta(CA,X,U,P,theta,theta_a,H_t,params)
    
    % Load parameters needed to create coefficients
    cA0 = params.cA0;
    stoich = params.stoich;
    THETA = params.THETA;
    p0 = params.p0;
    Pe_T = params.Pe_T;
    Da = params.Da;
    T0 = params.T0;
    Ta = params.Ta;
    K = params.K;
    beta = params.beta;
    Gamma = params.Gamma;
    N = params.N;
    a1 = params.a1;
    a2 = params.a2;
    a3 = params.a3;
    dt = params.dt;
    dz = params.dz;
    nodes = params.nodes;

    % Calculate enthalpy coefficient such that H*theta = total enthalpy
    AA = THETA'*a1 + THETA'*a2*T0*theta + THETA'*a3*T0^2*theta.^2  ...
            - (stoich'*a1+(stoich'*a2*T0)*theta+(stoich'*a3*T0^2) ...
                *theta.^2).*(X/stoich(1));
    H = (cA0/K)*AA./U;
% +++++ REVISION (11/29/25) +++++++++++++++++++++++++++++++++++++++++
    % Flux limiter
    UH = H.*theta;      % Calculate transport variable
    r = zeros(nodes,1); % Flux limiter variable
    r(2:nodes-1) = (UH(2:nodes-1)-UH(1:nodes-2)) ...
                    ./(UH(3:nodes)-UH(2:nodes-1));
    phi = max(zeros(nodes,1),min(ones(nodes,1),r));
    phi_1 = 1 - phi; 

    % Compute stencil coefficients
    W = phi.*-U.*H/(2*dz) - (1/Pe_T)*(1/dz^2) - phi_1.*(U.*H)/dz;
    C = H/dt + (1/Pe_T)*(2/dz^2) + Da*exp(beta./theta).*(CA) ...
        .*(cA0*(stoich'*a1+(stoich'*a2*T0/2)...
        *theta+((1/3)*stoich'*a3*T0^2)*theta.^2)/K)./1 + N ...
        + phi_1.*(U.*H)/dz;
    E = phi.*U.*H/(2*dz) - (1/Pe_T)*(1/dz^2);
% +++++ REVISION (11/29/25) +++++++++++++++++++++++++++++++++++++++++
    % Return coefficients in convenient form for using spdiags
    diag = [W,C,E];

    % Generate A at internal nodes
    A = spdiags(diag,[-1,0,1],nodes,nodes+1);
    A = A(:,1:nodes);

    % Clear boundary rows
    A(1,:) = 0; %A(nodes,:) = 0;

    % Left flux boundary
    A(1,1) = 0*(U(1)*H(1) - 2/(Pe_T*dz)) + 1;
    % Right flux boundary
    %A(nodes,nodes) = 1/dz;
    %A(nodes,nodes-1) = -1/dz;

    % Forcing function with appropriate boundary conditions

    % Internal nodes
    f = H_t/dt - Da*exp(beta./theta).*(CA).*Gamma./1 + N*theta_a*Ta/T0 ...
            + 0*(p0/K/T0)*U.*([P(2:nodes);0]-[0;P(1:nodes-1)])/(2*dz);
    f(1) = 1 - 2/(Pe_T*dz);     % Left
    %f(nodes) = 0;               % Right

end
