function [lastBestGuess, lastBestCov, G_hat, IG_hat, VarG_hat, VarIG_hat] = applyPF(pf, time, noisyMeasure, meal, total_ins, varargin)

ip = inputParser;

addParameter(ip,'InitialCondition_x0', [] );
addParameter(ip,'InitialCondition_sigma0', [] );
addParameter(ip,'processNoiseVariance', []);
addParameter(ip,'measurementNoiseVariance', []);
addParameter(ip,'parameterStructure', []);
addParameter(ip,'predictionHorizon', 30);

% parse the input arguments
parse(ip,varargin{:});

x0 = ip.Results.InitialCondition_x0;
sigma0 = ip.Results.InitialCondition_sigma0;
sigma_u = sqrt(ip.Results.processNoiseVariance);
sigma_v = sqrt(ip.Results.measurementNoiseVariance);
mP = ip.Results.parameterStructure;
PH = ip.Results.predictionHorizon;
 

% initialize array
lastBestGuess = zeros(length(x0),length(noisyMeasure));
lastBestGuess(8,1) = noisyMeasure(1);
lastBestGuess(9,1) = noisyMeasure(1);
lastBestCov = zeros(length(x0),length(noisyMeasure));

IG_hat = zeros(length(noisyMeasure),PH);
G_hat = zeros(length(noisyMeasure),PH);
VarIG_hat = zeros(length(noisyMeasure),PH);
VarG_hat = zeros(length(noisyMeasure),PH);

stateCorrected = zeros(1,length(x0));
covCorrected = zeros(length(x0),length(x0));

%% real-time filtering

for k = 1:length(time)
    
    %-- prediction step --%
   
    [statePred, covPred] = predict(pf, meal(k), total_ins(k), hms(time(k)), sigma_u, mP); 
    
    if mod(k-1,5/mP.TS) == 0  
        indexMeasure = (k+(5/mP.TS)-1)/(5/mP.TS);
        CGM_measure = noisyMeasure(indexMeasure);
    else
        CGM_measure = NaN;
    end
        
        
    %-- correction step --%
    if ~isnan(CGM_measure)
        [stateCorrected, covCorrected] = correct(pf, CGM_measure, sigma_v);
    else
        stateCorrected = statePred;
        covCorrected = covPred;
    end
    
    if mod(k-1,5/mP.TS) == 0
        lastBestGuess(:,indexMeasure) = stateCorrected(1:length(x0));
        lastBestCov(:,indexMeasure) = diag(covCorrected);
    end
     
    %-- k-step ahead prediction --%   
    if (k+PH <= length(time)) && (mod(k-1,5/mP.TS) == 0)
        
        pfPred = copy(pf);
        mPpred = mP;

        for p = 1:PH

            [stepAheadPrediction, covAheadPrediction] = predict(pfPred, meal(k+p), total_ins(k+p), hms(time(k+p)), sigma_u, mPpred);
            G_hat(indexMeasure,p) = stepAheadPrediction(8);
            IG_hat(indexMeasure,p) = stepAheadPrediction(9);
            
            PredVariance = diag(covAheadPrediction);
            VarG_hat(indexMeasure,p) = PredVariance(8);
            VarIG_hat(indexMeasure,p) = PredVariance(9);
             
        end
        
    end

    
end

end
