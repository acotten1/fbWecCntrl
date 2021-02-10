%% Clear, close all, clc

clear all; clc;

%% Step 0.0) Set up parameters
folder = fileparts(which('waveBotComp'));

preserveVolumesStudy = 0;

rho = 1025;     % water density

% Set scale
scale = 1; % R2 = 0.88m --> scale = 1
           % R2 = 2.20m --> scale = 2.2

% WaveBot dimensions (in metres)
r1 = 0.35*scale; % These are as per diagram
r2 = 0.88*scale;
t1 = 0.16*scale;
t2 = 0.53*scale;
t3 = 0.2*scale; % height of the cylindrical part above water

% Compute WaveBot volume for comparison with the compressible volume set for the waveBotComp.
vol_cyl = pi*(r2^2*(t1+t3)); % Volume of total cylindrical part (both above and below water)
vol_cyl_aboveWL = pi*(r2^2*(t3)); vol_cyl_belowWL = pi*(r2^2*(t1));
vol_truncCone = pi*((1/3)*(r2^2 + r1*r2 + r1^2)*(t2-t1));
vol_total = vol_cyl + vol_truncCone;
vol_aboveWL = vol_cyl_aboveWL; vol_belowWL = vol_cyl_belowWL + vol_truncCone;

if preserveVolumesStudy
    r1 = 0.35; t1 = 0.16; % r1 and t1 are to be left constant.
    
    t2 = 0.53; % Choose new value of t2 - i.e. the z-location of the CDOF
    % Define coefficients of quadratic equation to solve for new r2 that
    % preserves volume below waterline.
    a = pi*(t1 + (t2-t1)/3);
    b = pi*r1*(t2-t1)/3;
    c = pi*r1^2*(t2-t1)/3 - vol_belowWL;
    r2 = max(roots([a b c])); % One soln will be +ve, the other -ve.
    % Compute new t3 that preserves volume above waterline (though this is not actually used for the mass calculation at present).
    t3 = vol_aboveWL/(pi*r2^2);
    
end

% PTO damping value
dpto = FroudeScale.DampingRotational(0.1,2*10^5);

%% 0.2) Create frequency domain compuations for WaveBot with and without CDOF
h = 5;         % water depth
T = FroudeScale.Time(0.1,2:0.1:16);    % wave periods
T = 1./[0.25:0.005:0.8];

% Case 1) - Rigid body
fdComp_Rigid = waveBotComp(rho, h, T, r1, r2, t1, t2, t3, 'motions', [0 0 1 0 0 0]);

% Look at the WaveBot geometry
figure;
plot(fdComp_Rigid.Bodies.PanelGeo);
axis equal;
grid on;
box on;
set(gca, 'view', [10 20]);

% Case 2) - WaveBot with a Compressible Degree Of Freedom (CDOF) defined over its base
fdComp_CDOF = waveBotComp(rho, h, T, r1, r2, t1, t2, t3, 'motions', [0 0 1 0 0 0 1], 'cdofVol',vol_total*25);

% Plot heave motions
motions_Rigid = fdComp_Rigid.Motions;
motions_CDOF = fdComp_CDOF.Motions;
figure; hold on;
plot(1./T,abs(squeeze(motions_Rigid(:,1,1))))
plot(1./T,abs(squeeze(motions_CDOF(:,1,1))))
plot(1./T,abs(squeeze(motions_CDOF(:,1,2))))
legend('Heave, Rigid','Heave, CDOF','GenMode, CDOF')
xlabel('Frequency (Hz)')
ylabel('RAO Magnitude')

% Plot power curves               
Dpto = dpto;      % set the 1 x 1 PTO matrix
fdComp_Rigid.SetDpto(Dpto);  % set it on the FD comp.
Dpto = zeros(2,2); % set the 2 x 2 PTO matrix
Dpto(1,1) = dpto; % Power is extracted through waveBot's heave motion
fdComp_CDOF.SetDpto(Dpto);  % set it on the FD comp.

figure; hold on;
plot(1./T, fdComp_Rigid.PowerRAO);
plot(1./T, fdComp_CDOF.PowerRAO);
grid on;
legend('Rigid','with CDOF')
xlabel('Frequency (Hz)');
ylabel('Power RAO (kW/m^2)')

save('C:\Users\AlfredCotten\Desktop\MATLAB\snl_fbWecCntrlFork_CDOF\CDOF\data\waveBot_r0.88_CDOF.mat','fdComp_Rigid','fdComp_CDOF')

