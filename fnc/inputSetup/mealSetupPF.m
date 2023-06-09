function [meal, mealDelayed] = mealSetupPF(data,model,modelParameters,environment)
% function  mealSetup(data,model,modelParameters)
% Generates the vector containing the CHO intake events to be used to
% simulate the physiological model.
%
% Inputs:
%   - data: a timetable which contains the data to be used by the tool;
%   - model: a structure that contains general parameters of the
%   physiological model;
%   - modelParameters: a struct containing the model parameters;
%   - environment: a structure that contains general parameters to be used
%   by ReplayBG.
% Outputs:
%   - meal: is a vector containing the carbohydrate intake at each time
%   step [mg/min*kg];
%   - mealDelayed: is a vector containing the carbohydrate intake at each time
%   step delayed by beta min [mg/min*kg].
%

            
    %Initialize the meal vector
    meal = zeros(model.TIDSTEPS,1);


    %Set the meal vector
    for time = 1:length(0:model.YTS:(model.TID-1))
        meal((1+(time-1)*(model.YTS/model.TS)):(time*(model.YTS/model.TS))) = data.CHO(time)*1000/modelParameters.BW; %mg/(kg*min)
    end


    %Add delay of main meal absorption
    mealDelay = round(modelParameters.beta/model.TS);
    mealDelay = round(mealDelay/5)*5;
    mealDelayed = [zeros(mealDelay,1); meal];
    mealDelayed = mealDelayed(1:model.TIDSTEPS);
    meal = meal(1:model.TIDSTEPS);

            
    end
    
    
end