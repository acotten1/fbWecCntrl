%% Clear, close all, clc

clear all; clc;

%% 0) Set up parameters
folder = fileparts(which('waveBotComp'));

rho = 1025;     % water density
h = 1000;         % water depth (Select large number to approximate infinite depth)
T = 2*pi./[0.25:0.01:2]; % Wave periods
pAtm = 101325; % Atmospheric pressure(in Pascals)
adInd = 1.4; % Adiabatic index for air
g = 9.80665; % gravitational acceleration

% cylinder dimensions (from Kurniawan paper)
radius = 5;
draft = 10;
height = 20;

%% 1) Replicate cylinder cases from SAND report (and Kurniawan paper)
% Make frequency-domain computations with three cases for comparison with SAND report:
% 1) rigid, 2) CDOF with V10 = 1000m^3 and 3) CDOF with V10 = 1500m^3
% Case 1) - Rigid body
fdComp_Rigid = cylinderComp(rho, h, T, radius, draft, height, 'motions',[0 0 1 0 0 0],'panelSize',radius/3);

% Cases 2) and 3) - Cylinder with Compressible Degree Of Freedom (CDOF)
% Set CDOF volumes and compute the additional stiffness terms required:
K77_1000 = (pi*radius^2)^2*adInd*(pAtm + rho*g*draft)/1000;
K77_1500 = (pi*radius^2)^2*adInd*(pAtm + rho*g*draft)/1500;
fdComp_CDOF1 = cylinderComp(rho, h, T, radius, draft, height, 'motions',[0 0 1 0 0 0 1],'panelSize',radius/3,'cdofStiffness',K77_1000);
fdComp_CDOF2 = cylinderComp(rho, h, T, radius, draft, height, 'motions',[0 0 1 0 0 0 1],'panelSize',radius/3,'cdofStiffness',K77_1500);

% view the cylinder
figure;
plot(fdComp_CDOF1.Bodies.PanelGeo);
axis equal;
grid on;
box on;
set(gca, 'view', [10 20]);

% Plot heave motions
motions_Rigid = fdComp_Rigid.Motions;
motions_CDOF1 = fdComp_CDOF1.Motions;
motions_CDOF2 = fdComp_CDOF2.Motions;
figure; hold on;
plot(2*pi./T,abs(squeeze(motions_Rigid(:,1,1))))
plot(2*pi./T,abs(squeeze(motions_CDOF1(:,1,1))))
plot(2*pi./T,abs(squeeze(motions_CDOF2(:,1,1))))
legend('Rigid','CDOF with V10 = 1000m^3','CDOF with V10 = 1500m^3')
xlabel('Frequency (rad/s)')
ylabel('Heave RAO Magnitude')

