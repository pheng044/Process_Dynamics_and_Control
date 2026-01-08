% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025
% Calculate density field for non-dimensionalized PFR equation
% -------------------------------------------------------------------

function rho = density(X,theta,P,params)
    
    % Load parameters needed to calculate density
    cA0 = params.cA0;
    stoich = params.stoich;
    THETA = params.THETA;
    epsilon = params.epsilon;
    MW = params.MW;
    rho_0 = params.rho_0;

    
    % Calculate, return
    rho = (cA0./(1+epsilon*X)).*(sum(THETA.*MW) ...
            - sum(stoich.*MW)*X/stoich(1)).*P./theta/rho_0;

end
