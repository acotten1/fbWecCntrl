%% Clear, close all, clc

clear all; clc;

%% Step 0.0) Set up parameters
folder = fileparts(which('waveBotComp'));

rho = 1025;     % water density

% Set scale
scale = 1; % R2 = 0.88m --> scale = 1
           % R2 = 2.20m --> scale = 2.2
           
%%% Vary scale %%%
% scaleRange = [1 2 5 10 20 35 50 100];

% WaveBot dimensions (in metres)  
r2 = 0.88*scale;
t1 = 0.16*scale;
t2 = 0.53*scale;
t3 = 0.2*scale; % height of the cylindrical part above water

%%% Vary R1 %%%
% r1LB = 0.1; r1UB = r2;
% numVals = 12;
% r1range = r1LB:(r1UB-r1LB)/(numVals-1):r1UB;

%%% Vary T2 %%%
% t2LB = 0.53; t2UB = 4;
% numVals = 12;
% t2range = t2LB:(t2UB-t2LB)/(numVals-1):t2UB;

%%% Vary CDOF volume %%%
r1 = sqrt(r2^2*1.0);
cdofVolRangeMult = [1 2 5 10 20 35 50];
           
for i = 1:length(cdofVolRangeMult)

% WaveBot dimensions (in metres)
% r1 = r1range(i)*scale; % These are as per diagram
% t2 = t2range(i)*scale;

% Compute WaveBot volume for comparison with the compressible volume set for the waveBotComp.
vol_cyl = pi*(r2^2*(t1+t3)); % Volume of total cylindrical part (both above and below water)
vol_cyl_aboveWL = pi*(r2^2*(t3)); vol_cyl_belowWL = pi*(r2^2*(t1));
vol_truncCone = pi*((1/3)*(r2^2 + r1*r2 + r1^2)*(t2-t1));
vol_total(i) = vol_cyl + vol_truncCone;
vol_aboveWL = vol_cyl_aboveWL; vol_belowWL(i) = vol_cyl_belowWL + vol_truncCone;

%% 0.2) Create frequency domain compuations for WaveBot with and without CDOF
h = 5*scale;         % water depth
T = FroudeScale.Time(scale,1./[0.25:0.005:0.8]);

f = 1./T';

pAtm = 101325; % Atmospheric pressure (in Pascals)
g = 9.80665; % gravitational acceleration
p0 = (pAtm + rho*g*t2);
adInd = 1.4; % Adiabatic index for air
S1 = pi*r1^2; S2 = pi*r2^2;

cdofVol = cdofVolRangeMult(i)*vol_belowWL(i);

% Case 2) - WaveBot with a Compressible Degree Of Freedom (CDOF) defined over its base
fdComp_CDOF = waveBotComp(rho, h, T, r1, r2, t1, t2, t3, 'motions', [0 0 1 0 0 0 1], 'cdofVol',cdofVol);

% %% Compute lower limit on natural frequency (V0 --> Inf)
% f = 1./T';
% pAtm = 101325; % Atmospheric pressure (in Pascals)
% g = 9.80665; % gravitational acceleration
% p0 = (pAtm + rho*g*t2);
% adInd = 1.4; % Adiabatic index for air
% S1 = pi*r1^2; S2 = pi*r2^2;
% 
% a = (fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*fdComp_CDOF.A(:,2,2) - fdComp_CDOF.A(:,1,2).^2;
% b = -((fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*(rho*g*S1) + rho*g*S2*fdComp_CDOF.A(:,2,2) + fdComp_CDOF.B(:,1,1).*fdComp_CDOF.B(:,2,2) ...
%     - 2*rho*g*S1*fdComp_CDOF.A(:,1,2) - fdComp_CDOF.B(:,1,2).^2);
% c = rho*g*S1*(rho*g*S2 - rho*g*S1);
% wresSquared = (- b - sqrt(b.^2 - 4*a*c))./(2*a); % Solve for omega squared (smallest of the two solutions is of interest here)
% wres = sqrt(wresSquared); % Only positive solution is of interest here
% fres = wres/(2*pi);
% [val ind] = min(abs(fres - f)); % 
% fres(ind)
% 
% %% Compute lower limit on natural frequency using simpler formula
% wres = sqrt(rho*g*(S2-S1)./(fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)-fdComp_CDOF.A(:,1,2)));
% fres = wres/(2*pi);
% [val ind] = min(abs(fres - f)); % 
% fres(ind)

%% Compute natural frequency (no limits taken)
a = (fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*fdComp_CDOF.A(:,2,2) - fdComp_CDOF.A(:,1,2).^2;
b = -((fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*(rho*g*S1 + S1^2*adInd*p0/cdofVol) + rho*g*S2*fdComp_CDOF.A(:,2,2) + fdComp_CDOF.B(:,1,1).*fdComp_CDOF.B(:,2,2) ...
    - 2*rho*g*S1*fdComp_CDOF.A(:,1,2) - fdComp_CDOF.B(:,1,2).^2);
c = rho*g*S1*(rho*g*S2 - rho*g*S1) + rho*g*S2*S1^2*adInd*p0/cdofVol;
wresSquared = (- b - sqrt(b.^2 - 4*a*c))./(2*a); % Solve for omega squared (smallest of the two solutions is of interest here)
wres = sqrt(wresSquared); % Only positive solution is of interest here
fres = wres/(2*pi);
[val ind] = min(abs(fres - f)); % 
natFreqComplic(i) = fres(ind)

% wresSquared2 = (- b - sqrt(b.^2 - 4*a*c)); % Solve for omega squared (smallest of the two solutions is of interest here)
% wres2 = sqrt(wresSquared2); % Only positive solution is of interest here
% fres2 = wres2/(2*pi);
% [val ind] = min(abs(fres2 - f)); % 
% natFreqComplic2(i) = fres2(ind)

%% Compute natural frequency with limit V0 --> Inf taken.
a = (fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*fdComp_CDOF.A(:,2,2) - fdComp_CDOF.A(:,1,2).^2;
b = -((fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*(rho*g*S1) + rho*g*S2*fdComp_CDOF.A(:,2,2) + fdComp_CDOF.B(:,1,1).*fdComp_CDOF.B(:,2,2) ...
    - 2*rho*g*S1*fdComp_CDOF.A(:,1,2) - fdComp_CDOF.B(:,1,2).^2);
c = rho*g*S1*(rho*g*S2 - rho*g*S1);
wresSquared = (- b - sqrt(b.^2 - 4*a*c))./(2*a); % Solve for omega squared (smallest of the two solutions is of interest here)
wres = sqrt(wresSquared); % Only positive solution is of interest here
fres = wres/(2*pi);
[val ind] = min(abs(fres - f)); % 
natFreqComplicV0Inf(i) = fres(ind)

% wresSquared2 = (- b - sqrt(b.^2 - 4*a*c)); % Solve for omega squared (smallest of the two solutions is of interest here)
% wres2 = sqrt(wresSquared2); % Only positive solution is of interest here
% fres2 = wres2/(2*pi);
% [val ind] = min(abs(fres2 - f)); % 
% natFreqComplicV0Inf2(i) = fres2(ind)

%% Compute natural frequency using simpler formula
r = -rho*g/(S1^2*adInd*p0/cdofVol + rho*g);
wres = sqrt(rho*g*(S2+S1*r)./(fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)+fdComp_CDOF.A(:,1,2)*r));
fres = wres/(2*pi);
[val ind] = min(abs(fres - f)); % 
natFreqSimple(i) = fres(ind);

%% Rigid body resonant frequency
wresRB = sqrt(rho*g*S2./(fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)));
fresRB = wresRB/(2*pi);
[val ind] = min(abs(fresRB - f)); % 
natFreqRB(i) = fresRB(ind)

end

%% Plot resonant frequency against S1 value
% S1range = pi*r1range.^2; 
% figure; hold on; plot(S1range/S2,natFreqRB); plot(S1range/S2,natFreqComplic); plot(S1range/S2,natFreqComplicV0Inf);
% xlabel('S_1/S_2');ylabel('Natural frequency (Hz)');legend('Rigid Body','CDOF, V_0 = submerged volume','CDOF, V_0 --> Inf')
 
% figure; hold on; plot(t2range,natFreqRB); plot(t2range,natFreqComplic); plot(t2range,natFreqComplicV0Inf);
% xlabel('T_2 (m)');ylabel('Natural frequency (Hz)');legend('Rigid Body','CDOF, V_0 = submerged volume','CDOF, V_0 --> Inf')

figure; hold on; plot(cdofVolRangeMult,natFreqRB); plot(cdofVolRangeMult,natFreqComplic); plot(cdofVolRangeMult,natFreqSimple);
xlabel('CDOF vol as multiplier of submerged WaveBot volume');ylabel('Natural frequency (Hz)');legend('Rigid Body','with CDOF, full equation','with CDOF, simplified equation')

figure; hold on; plot(cdofVolRangeMult,natFreqRB-natFreqComplic); plot(cdofVolRangeMult,natFreqRB-natFreqSimple);
xlabel('CDOF vol as multiplier of submerged WaveBot volume');ylabel('Decrease in natural frequency due to CDOF (Hz)');legend('Full equation','Simplified equation')

% figure; hold on; plot(scaleRange,natFreqRB); plot(scaleRange,natFreqComplic);
% xlabel('scale compared to original WaveBot');ylabel('Natural frequency (Hz)');legend('Rigid Body','with CDOF')

% figure; hold on; plot(S1range/S2,natFreqComplic); plot(S1range/S2,natFreqSimple); plot(S1range/S2,natFreqComplicV0Inf);
% xlabel('S_1/S_2');ylabel('Natural frequency (Hz)');legend('Full equation','Simplified equation','Full equation - V0 --> Inf')
% 
% figure; hold on; plot(S1range/S2,natFreqComplicV0Inf2);
% xlabel('S_1/S_2');ylabel('Natural frequency (Hz)');legend('Full equation - V0 --> Inf, numerator only')




