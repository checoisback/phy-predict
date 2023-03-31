function likelihood = giPFMeasurementLikelihoodFcn(particles, measurement, sigma_v)

% The measurement is nineth state. Get all measurement hypotheses from particles
predictedMeasurement = particles(9,:);

% Assume the ratio of the error between predicted and actual measurements
% follow a Gaussian distribution with zero mean, variance sigma
mu = 0; % mean
sigma = sigma_v*eye(1);
likelihood = normpdf(predictedMeasurement, measurement, sigma);

end