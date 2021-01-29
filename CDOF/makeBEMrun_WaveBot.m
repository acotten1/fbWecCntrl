%% Clear, close all, clc

clear all; clc;

%% Step 0.0) Set up parameters
folder = fileparts(which('waveBotComp'));

rho = 1025;     % water density

% Set scale
scale = 2.2; % R2 = 0.88m --> scale = 1
           % R2 = 2.20m --> scale = 2.2

% WaveBot dimensions (in metres)
r1 = 0.35*scale; % These are as per diagram
r2 = 0.88*scale;
t1 = 0.16*scale;
t2 = 0.53*scale;
t3 = 0.2*scale; % height of the cylindrical part above water

% PTO damping value
dpto = FroudeScale.DampingRotational(0.1,2*10^5);

%% 0.2) Create frequency domain compuations for WaveBot with and without CDOF
h = 5;         % water depth
T = FroudeScale.Time(0.1,2:0.2:16);    % wave periods
T = [33.3333333333333;30.3030303030303;27.5482093663912;25.0626566416040;22.8310502283105;20.7468879668050;18.9035916824197;17.1821305841924;15.6250000000000;14.2247510668563;12.9366106080207;11.7647058823529;10.7066381156317;9.73709834469328;8.85739592559787;8.05801772763900;7.32600732600733;6.66666666666667;6.06428138265616;5.51571980143409;5.01756146512795;4.56412596987677;4.15282392026578;3.77786173026067;3.43642611683849;3.12597686777118;2.84333238555587;2.58665287118469;2.35294117647059;2.14041095890411;1.94704049844237;1.77116542685087;1.61134386077989;1.46563095412575;1.33333333333333];

% Case 1) - Rigid body
fdComp_Rigid = waveBotComp(rho, h, T, r1, r2, t1, t2, t3, 'motions', [0 0 1 0 0 0]);

% Look at the WaveBot geometry
figure;
plot(fdComp_Rigid.Bodies.PanelGeo);
axis equal;
grid on;
box on;
set(gca, 'view', [10 20]);

% Compute WaveBot volume for comparison with the compressible volume set for the waveBotComp.
vol_cyl = pi*(r2^2*(t1+t3)); % Volume of total cylindrical part (both above and below water)
vol_truncCone = pi*((1/3)*(r2^2 + r1*r2 + r1^2)*(t2-t1));
vol_total = vol_cyl + vol_truncCone;

% Case 2) - WaveBot with a Compressible Degree Of Freedom (CDOF) defined over its base
fdComp_CDOF = waveBotComp(rho, h, T, r1, r2, t1, t2, t3, 'motions', [0 0 1 0 0 0 1], 'cdofVol',vol_total*15);

% Plot heave motions
motions_Rigid = fdComp_Rigid.Motions;
motions_CDOF = fdComp_CDOF.Motions;
figure; hold on;
plot(2*pi./T,abs(squeeze(motions_Rigid(:,1,1))))
plot(2*pi./T,abs(squeeze(motions_CDOF(:,1,1))))
plot(2*pi./T,abs(squeeze(motions_CDOF(:,1,2))))
legend('Heave, Rigid','Heave, CDOF','GenMode, CDOF')
xlabel('Frequency (rad/s)')
ylabel('RAO Magnitude')

% Plot power curves               
Dpto = dpto;      % set the 1 x 1 PTO matrix
fdComp_Rigid.SetDpto(Dpto);  % set it on the FD comp.
Dpto = zeros(2,2); % set the 2 x 2 PTO matrix
Dpto(1,1) = dpto; % Power is extracted through waveBot's heave motion
fdComp_CDOF.SetDpto(Dpto);  % set it on the FD comp.

figure; hold on;
plot(T, fdComp_Rigid.PowerRAO);
plot(T, fdComp_CDOF.PowerRAO);
grid on;
legend('Rigid','with CDOF')
xlabel('Period (s)');
ylabel('Power RAO (kW/m^2)')

save('C:\Users\AlfredCotten\Desktop\MATLAB\snl_fbWecCntrlFork_CDOF\CDOF\data\waveBot_r2.20_CDOF.mat','fdComp_Rigid','fdComp_CDOF')

