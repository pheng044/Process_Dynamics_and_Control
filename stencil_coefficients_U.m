% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025 (Steady State Solver)
% 14 Nov 2025 - 6 Dec 2025 (Transient Solver)
% 1D stencil coefficients for non-dimensionalized PFR problem
%  -> momentum equation
% -------------------------------------------------------------------

function [A,f] = stencil_coefficients_U(U,rho,P,U_mom_t,SS_U,params)
    
    % Load parameters needed to create coefficients
    Re = params.Re;
    Eu = params.Eu;
    dt = params.dt;
    dz = params.dz;
    fric = params.fric;
    nodes = params.nodes;
 

% +++++ REVISION (11/29/25) +++++++++++++++++++++++++++++++++++++++++
    % If momentum solve is set to steady state, remove the time 
    % derivative term
    if SS_U == true
        dt = inf;
    end
    
    % Flux limiter calculation
    U_mom = rho.*U;     % Calculate transport variable 
    %r = zeros(nodes,1); % Flux limiter variable
    
    %r(2:nodes-1) = (U_mom(2:nodes-1)-U_mom(1:nodes-2)) ... 
    %                ./(U_mom(3:nodes)-U_mom(2:nodes-1));
    %phi = max(zeros(nodes,1),min(ones(nodes,1),r));
    phi = 0*ones(nodes,1);
    phi_1 = 1 - phi; 

    % Compute stencil coefficients
    W = phi.*-U_mom/(2*dz) - (1/Re)/dz^2 - phi_1.*(U_mom)/dz;
    C = rho/dt + (1/Re)*2/dz^2 + phi_1.*(U_mom)/dz ... 
            + 0.5*Eu*fric.*U_mom;
    E = phi.*U_mom/(2*dz) - (1/Re)/dz^2;
% +++++ REVISION (11/29/25) +++++++++++++++++++++++++++++++++++++++++    
    % Return coefficients in convenient form for using spdiags
    diag = [W,C,E];

    % Generate A at internal nodes
    A = spdiags(diag,[-1,0,1],nodes,nodes+1);
    A = A(:,1:nodes);

    % Clear boundary rows
    A(1,:) = 0; A(nodes,:) = 0;

    % Left flux boundary
    A(1,1) = 1 + 0*(rho(1)*U(1) - 2/(Re*dz));
    % Right flux boundary
    A(nodes,nodes) = 1/dz;
    A(nodes,nodes-1) = -1/dz;
    
    %dpdz = Eu*([0;P(3:nodes);0]-[0;P(1:nodes-1)])/(2*dz);
    %dpdz = Eu*([0;P(2:nodes)]-[0;P(1:nodes-1)])/(dz);

    % Forcing function with appropriate boundary conditions
    f = U_mom_t/dt; %Internal
    f(1) = 1 - 0*2/(Re*dz);           % Left
    f(nodes) = 0;                   % Right

end