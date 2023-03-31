function particles = giParticleFilterStateFcn(particles, CHO, INS, time, sigma_u, mP) 

[numberOfStates, numberOfParticles] = size(particles);
    
% Time-propagate each particle

% Euler integration of continuous-time dynamics x'=f(x) with sample time dt
dt = mP.TS; % Sample time
for kk = 1:numberOfParticles
    particles(:,kk) = particles(:,kk) + giStateFcnContinuous(particles(:,kk), CHO, INS, time, mP)*dt;
end

% Add Gaussian noise with variance processNoise
processNoise = sigma_u'.*eye(numberOfStates);
particles = particles + processNoise * randn(size(particles));

end

function dxdt = giStateFcnContinuous(x, CHO, INS, time, mP)
% giStateFcnContinuous Evaluate the glucose insulin system
% Qsto1ss = 0
Qsto1 = x(1); 

% Qsto2ss = 0
Qsto2 = x(2);

% Qgutss = 0
Qgut = x(3); 

% Isc1ss = u2ss / ( ka1 + kd )
Isc1 = x(4); % mU/kg

% Isc2ss = kd / ka2 * u2ss / ( ka1 + kd )
Isc2 = x(5); % mU/kg

% Ipss = ka1 / ke * u2ss / ( ka1 + kd ) + ka2 / ke * kd / ka2 * u2ss / ( ka1 + kd )
Ip = x(6); % mU/kg

% Compute the basal plasmatic insulin
Ipb = (mP.ka1/mP.ke)*(mP.u2ss)/(mP.ka1+mP.kd) + (mP.ka2/mP.ke)*(mP.kd/mP.ka2)*(mP.u2ss)/(mP.ka1+mP.kd); %from eq. 5 steady-state    

% Xss = Xb = 0
X = x(7); 

% Gss = Gb
G = x(8);

% IGss = IGb
IG = x(9); 


% Compute the mP state at time k using backward Euler method

if(time<4 || time >= 17)
    SI = mP.SID;
else
    if(time >=4 && time < 11)
        SI = mP.SIB;
    else
        SI = mP.SIL;
    end
end

%Compute the hypoglycemic risk
risk = computeHypoglycemicRisk(G,mP);
rhoRisk = (1+risk);

dxdt(1,1) = -mP.kgri*Qsto1+CHO;
dxdt(2,1) = mP.kgri*Qsto1-mP.kempt*Qsto2;
dxdt(3,1) = mP.kempt*Qsto2-mP.kabs0*Qgut;

Ra = mP.f*mP.kabs0*Qgut;

dxdt(4,1) = -mP.kd*Isc1+INS;
dxdt(5,1) = mP.kd*Isc1-mP.ka2*Isc2;
dxdt(6,1) = mP.ka2*Isc2-mP.ke*Ip;

dxdt(7,1) = -mP.p2*(X-(SI/mP.VI)*(Ip-Ipb));
dxdt(8,1) = -((mP.SG+rhoRisk*X)*G)+mP.SG*mP.Gb+Ra/mP.VG;
dxdt(9,1) = -(1/mP.alpha)*(IG-G);

end