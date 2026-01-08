% -------------------------------------------------------------------
% Patrick Heng
% 15 Feb 2025 - 6 Apr 2025
% Calculate pressure field for non-dimensionalized PFR equation
% -------------------------------------------------------------------

function P = pressure(X,U,theta,params)

    epsilon = params.epsilon;
    p_max = params.p_max;
    P = (1 + epsilon*X).*theta./U;
    P(P>p_max) = p_max;

end