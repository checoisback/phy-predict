function risk = computeHypoglycemicRisk(G,mP)
% function  computeHypoglycemicRisk(G,mP)
% Computes the hypoglycemic risk as in Visentin et al., JDST, 2018.
%
% Inputs:
%   - G: the glucose concentration;
%   - mP: is a struct containing the model parameters.
% Output:

    %Setting the risk model threshold
    Gth = 60;
    
    %Compute the risk
    Gb = mP.Gb;
    risk = 10*(log(G)^mP.r2 - log(Gb)^mP.r2)^2*(G<Gb & G>=Gth) + ...
        10*(log(Gth)^mP.r2 - log(Gb)^mP.r2)^2*(G<Gth);
    risk = abs(risk);