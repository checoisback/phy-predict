function predictionResults = getPhyPrediction(identifiedModelParameters,testData,PH)
% function getPhyPrediction.m
% it computes the glucose predicted profiles by using a
% non-linear physiological model of glucose-insulin dynamics integrated
% into a particle filter. 
%
% Inputs: 
% - identifiedModelParameters: a structure containing the identified model
% parameters described in https://github.com/gcappon/replay-bg.
% - testData: a timetable containing the following fields: Time, glucose,
% CHO, bolus_insulin and basal_insulin. 
% Of note, the measurements units must be in line with the ones described 
% in: https://github.com/gcappon/replay-bg.
% - PH: the prediction horizon (in minutes).
% Output:
% - predictionResults: a structure containing the test data 
% (predictionResults.data) and the predicted profile
% (predictionResults.ighat);

addpath(genpath('fnc'))


    %% set parameteres 
    Ts = 5; % minutes
    nParticles = 5000; % number of particles
    sigma_v = 25; % measurements noise variance
    sigma_u0 = [10 10 10 0.6 0.6 0.6 1e-6 10 10]; % process noise variance 

    mP = identifiedModelParameters; % identified model parameters
    data = testData; % data to predict ahead in time

    % set physiological model environment
    model.TS = 1; % physiological model sampling time. [minute]
    model.YTS = 5; % raw data sampling time. [minute]
    model.TID = minutes(data.Time(end)-data.Time(1))+model.TS;  % from 1 to TID identify the model parameters. [min]
    model.TIDSTEPS = model.TID/model.TS;    % from 1 to TID identify the model parameters. [integration steps]
    model.TIDYSTEPS = model.TID/model.YTS;  % total identification simulation time [sample steps]
    mP.TS = model.TS;

    % prepare input data
    [bolus, basal, bolusDelayed, basalDelayed] = insulinSetupPF(data, model, mP);
    [cho, choDelayed] = mealSetupPF(data, model, mP);
    meal = choDelayed;
    total_ins = (bolusDelayed + basalDelayed);
    timeData = data.Time(1):minutes(1):(data.Time(1) + minutes(length(cho)-1));


    % initial state
    x0 = [ 0; ...% Qsto1(0)
           0; ...% Qsto2(0)
           0; ...% Qgut(0)
           mP.u2ss/(mP.ka1+mP.kd); ... % Isc1(0)
           (mP.kd/mP.ka2)*mP.u2ss/(mP.ka1+mP.kd); ... % Isc2(0)
           (mP.ka1/mP.ke)*mP.u2ss/(mP.ka1+mP.kd) + (mP.ka2/mP.ke)*(mP.kd/mP.ka2)*mP.u2ss/(mP.ka1+mP.kd); ... % Ip(0)
           mP.Xpb; ... % X(0)
           data.glucose(1); ... % G(0)
           data.glucose(1); ... % IG(0)
           ];

    sigma0 = eye(length(x0));


    %% particle filter MATLAB implementation
    pf_v0 = particleFilter;
    pf_v0.StateEstimationMethod = 'mean'; 
    pf_v0.ResamplingMethod = 'systematic';
    pf_v0.StateTransitionFcn = @giParticleFilterStateFcn; % state update function
    pf_v0.MeasurementLikelihoodFcn = @giPFMeasurementLikelihoodFcn; % likelihood function
    
    stateBound = ones(length(x0),2);
    stateBound(:,1) = x0 - 0.03*x0;
    stateBound(:,2) = x0 + 0.03*x0;

    initialize(pf_v0, nParticles, stateBound);


    [lastBestGuess, lastBestCov, G_hat, IG_hat, VarG_hat, VarIG_hat] = applyPF(pf_v0, timeData, data.glucose, meal, total_ins, ...
        'InitialCondition_x0', x0,'InitialCondition_sigma0', sigma0,...
        'processNoiseVariance', sigma_u0,...
        'measurementNoiseVariance', sigma_v,...
        'parameterStructure', mP,...
        'predictionHorizon', PH);
    
    predicted_CGM = IG_hat(1:end-(PH/Ts),end);

    ighat = timetable(predicted_CGM,'RowTimes',data.Time(1+(PH/Ts):end));
    ighat.Properties.VariableNames = {'glucose'};
    
    
    predictionResults.dataTest = data;
    predictionResults.dataHat = ighat;
    
    
end