% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025 (Steady State Solver)
% 14 Nov 2025 - 6 Dec 2025 (Transient Solver)
% 1D stencil coefficients for non-dimensionalized PFR problem
%  -> concentration equation
% -------------------------------------------------------------------

function [A,f] = stencil_coefficients_CA(CA,U,theta,CA_t,in,params)

    % Load parameters needed to create coefficients
    Pe_M = params.Pe_M;
    Da = params.Da;
    beta = params.beta;
    cA0 = params.cA0;
    dt = params.dt;
    dz = params.dz;
    nodes = params.nodes;
% +++++ REVISION (11/29/25) +++++++++++++++++++++++++++++++++++++++++
    % Flux limiter 
    r = zeros(nodes,1);
    r(2:nodes-1) = (CA(2:nodes-1)-CA(1:nodes-2))...
                    ./(CA(3:nodes)-CA(2:nodes-1));
    phi = max(zeros(nodes,1),min(ones(nodes,1),r));
    phi_1 = 1 - phi; 

    % Compute stencil coefficients
    W = - phi.*U/(2*dz) - (1/Pe_M)/dz^2 - phi_1.*U/dz;
    C =  1/dt + (1/Pe_M)*2/dz^2 + Da*exp(beta./theta) + phi_1.*U/dz;
    E = phi.*U/(2*dz) - (1/Pe_M)/dz^2;
% +++++ REVISION (11/29/25) +++++++++++++++++++++++++++++++++++++++++

    % Return coefficients in convenient form for using spdiags
    diag = [W,C,E];
    
    % Generate A at internal nodes
    A = spdiags(diag,[-1,0,1],nodes,nodes+1);
    A = A(:,1:nodes);

    % Clear boundary rows
    A(1,:) = 0; 
    A(nodes,:) = 0;
    
    % Left flux boundary
    A(1,1) = U(1) - 2/(Pe_M*dz);
    % Right flux boundary
    A(nodes,nodes) = 1/dz;
    A(nodes,nodes-1) = -1/dz;
  
    % Forcing function with appropriate boundary conditions
    f = CA_t/dt;                % Internal
    f(1) = 1+in/cA0 - 2/(Pe_M*dz);     % Right
    f(nodes) = 0;               % Left
    
end
