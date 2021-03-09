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
h = 5*scale;         % water depth
T = FroudeScale.Time(0.1,2:0.1:16);    % wave periods
T = FroudeScale.Time(scale,1./[0.25:0.005:0.8]);

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
cdofVol = 25*vol_total;
fdComp_CDOF = waveBotComp(rho, h, T, r1, r2, t1, t2, t3, 'motions', [0 0 1 0 0 0 1], 'cdofVol',cdofVol);

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

%% Compute lower limit on natural frequency (V0 --> Inf)
f = 1./T';
pAtm = 101325; % Atmospheric pressure (in Pascals)
g = 9.80665; % gravitational acceleration
p0 = (pAtm + rho*g*t2);
adInd = 1.4; % Adiabatic index for air
S1 = pi*r1^2; S2 = pi*r2^2;

a = (fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*fdComp_CDOF.A(:,2,2) - fdComp_CDOF.A(:,1,2).^2;
b = -((fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*(rho*g*S1) + rho*g*S2*fdComp_CDOF.A(:,2,2) + fdComp_CDOF.B(:,1,1).*fdComp_CDOF.B(:,2,2) ...
    - 2*rho*g*S1*fdComp_CDOF.A(:,1,2) - fdComp_CDOF.B(:,1,2).^2);
c = rho*g*S1*(rho*g*S2 - rho*g*S1);
wresSquared = (- b - sqrt(b.^2 - 4*a*c))./(2*a); % Solve for omega squared (smallest of the two solutions is of interest here)
wres = sqrt(wresSquared); % Only positive solution is of interest here
fres = wres/(2*pi);
[val ind] = min(abs(fres - f)); % 
fres(ind)

%% Compute lower limit on natural frequency using simpler formula
wres = sqrt(rho*g*(S2-S1)./(fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)-fdComp_CDOF.A(:,1,2)));
fres = wres/(2*pi);
[val ind] = min(abs(fres - f)); % 
fres(ind)

%% Compute natural frequency (no limits taken)
a = (fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*fdComp_CDOF.A(:,2,2) - fdComp_CDOF.A(:,1,2).^2;
b = -((fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*(rho*g*S1 + S1^2*adInd*p0/cdofVol) + rho*g*S2*fdComp_CDOF.A(:,2,2) + fdComp_CDOF.B(:,1,1).*fdComp_CDOF.B(:,2,2) ...
    - 2*rho*g*S1*fdComp_CDOF.A(:,1,2) - fdComp_CDOF.B(:,1,2).^2);
c = rho*g*S1*(rho*g*S2 - rho*g*S1) + rho*g*S2*S1^2*adInd*p0/cdofVol;
wresSquared = (- b - sqrt(b.^2 - 4*a*c))./(2*a); % Solve for omega squared (smallest of the two solutions is of interest here)
wres = sqrt(wresSquared); % Only positive solution is of interest here
fres = wres/(2*pi);
[val ind] = min(abs(fres - f)); % 
fres(ind)

% % Make plots to sketch the dependence on each variable
% % S1
% S1range = 0:S2/100:S2; index = 50;
% a = (fdComp_CDOF.M(1,1)+fdComp_CDOF.A(index,1,1)).*fdComp_CDOF.A(index,2,2) - fdComp_CDOF.A(index,1,2).^2;
% b = -((fdComp_CDOF.M(1,1)+fdComp_CDOF.A(index,1,1)).*(rho*g*S1range + S1range.^2*adInd*p0/cdofVol) + rho*g*S2*fdComp_CDOF.A(index,2,2) + fdComp_CDOF.B(index,1,1).*fdComp_CDOF.B(index,2,2) ...
%     - 2*rho*g*S1range*fdComp_CDOF.A(index,1,2) - fdComp_CDOF.B(index,1,2).^2);
% c = rho*g*S1range.*(rho*g*S2 - rho*g*S1range) + rho*g*S2*S1range.^2*adInd*p0/cdofVol;
% wresSquared = (- b - sqrt(b.^2 - 4*a*c))./(2*a); % Solve for omega squared (smallest of the two solutions is of interest here)
% wres = sqrt(wresSquared); % Only positive solution is of interest here
% fres = wres/(2*pi);
% figure;plot(S1range/S2,fres);xlabel('S_1/S_2');ylabel('Natural frequency (Hz)')


%% Compute natural frequency using simpler formula
r = -rho*g/(S1^2*adInd*p0/cdofVol + rho*g);
wres = sqrt(rho*g*(S2+S1*r)./(fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)+fdComp_CDOF.A(:,1,2)*r));
fres = wres/(2*pi);
[val ind] = min(abs(fres - f)); % 
fres(ind)

% % Make plots to sketch the dependence on each variable
% % S1
% S1range = 0:S2/100:S2; index = 50;
% rrange = -rho*g/(S1^2*adInd*p0/cdofVol + rho*g);
% wres = sqrt(rho*g*(S2+S1range.*rrange)./(fdComp_CDOF.M(1,1)+fdComp_CDOF.A(index,1,1)+fdComp_CDOF.A(index,1,2)*rrange));
% fres = wres/(2*pi);
% figure;plot(S1range/S2,fres,'--k');xlabel('S_1/S_2');ylabel('Natural frequency (Hz)')

%% Rigid body resonant frequency
wresRB = sqrt(rho*g*S2./(fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)));
fresRB = wresRB/(2*pi);
[val ind] = min(abs(fresRB - f)); % 
fresRB(ind)



% %% 
% %%% S2 allowed to vary %%%
% x = 0.53; %T2
% A33 = 1500; 
% S2 = 0:0.02:10; 
% n=10; 
% yrb = sqrt(1025*9.81./(1025*x+A33./S2)); 
% ycdof = sqrt(1025*9.81*1.4*(101325+1025*9.81.*x)./((1025*x+A33./S2)*1.4.*(101325+1025*9.81.*x)+1025^2*9.81*n*x.^2));
% 
% figure;plot(S2,yrb,S2,ycdof);legend('RB','CDOF');xlabel('S2');ylabel('wres')
% figure;plot(S2,yrb-ycdof);xlabel('S2');ylabel('wresRB-wCDOF')
% 
% % Allow S2 to vary with full eqn
% S2 = 0:0.02:10; S1 = S2;
% t2 = 0.53;
% 
% for i = 1:length(S2)
% a = (fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*fdComp_CDOF.A(:,2,2) - fdComp_CDOF.A(:,1,2).^2;
% b = -((fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*(rho*g*S1(i) + S1(i)^2*adInd*(pAtm + rho*g*t2)/cdofVol) + rho*g*S2(i)*fdComp_CDOF.A(:,2,2) + fdComp_CDOF.B(:,1,1).*fdComp_CDOF.B(:,2,2) ...
%     - 2*rho*g*S1(i)*fdComp_CDOF.A(:,1,2) - fdComp_CDOF.B(:,1,2).^2);
% c = rho*g*S1(i)*(rho*g*S2(i) - rho*g*S1(i)) + rho*g*S2(i)*S1(i)^2*adInd.*(pAtm + rho*g*t2)/cdofVol;
% wresSquared = (- b - sqrt(b.^2 - 4*a*c))./(2*a); % Solve for omega squared (smallest of the two solutions is of interest here)
% wres = sqrt(wresSquared); % Only positive solution is of interest here
% fres = wres/(2*pi);
% [val ind] = min(abs(fres - f)); % 
% ycdofFull(i) = fres(ind);
% end
% 
% figure;plot(S2,yrb,S2,ycdofFull);legend('RB','CDOF');xlabel('S2');ylabel('wres')
% figure;plot(S2,yrb-ycdofFull);xlabel('S2');ylabel('wresRB-wCDOF')
% 
% %%% T2 allowed to vary %%%
% x = 0:0.02:100; %T2
% A33 = 1500;
% S2 = 2.4328;
% n = 10; 
% yrb = sqrt(1025*9.81./(1025*x+A33/S2)); 
% ycdof = sqrt(1025*9.81*1.4*(101325+1025*9.81.*x)./((1025*x+A33/S2)*1.4.*(101325+1025*9.81.*x)+1025^2*9.81*n*x.^2));
% 
% figure;plot(x,yrb,x,ycdof);legend('RB','CDOF');xlabel('T2');ylabel('wres')
% figure;plot(x,yrb-ycdof);xlabel('T2');ylabel('wresRB-wCDOF')
% 
% % Allow T2 to vary with full eqn
% t2 = 0:0.02:100;
% S2 = 2.4328; S1 = S2;
% 
% for i = 1:length(t2)
% a = (fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*fdComp_CDOF.A(:,2,2) - fdComp_CDOF.A(:,1,2).^2;
% b = -((fdComp_CDOF.M(1,1)+fdComp_CDOF.A(:,1,1)).*(rho*g*S1 + S1^2*adInd*(pAtm + rho*g*t2(i))/cdofVol) + rho*g*S2*fdComp_CDOF.A(:,2,2) + fdComp_CDOF.B(:,1,1).*fdComp_CDOF.B(:,2,2) ...
%     - 2*rho*g*S1*fdComp_CDOF.A(:,1,2) - fdComp_CDOF.B(:,1,2).^2);
% c = rho*g*S1*(rho*g*S2 - rho*g*S1) + rho*g*S2*S1^2*adInd.*(pAtm + rho*g*t2(i))/cdofVol;
% wresSquared = (- b - sqrt(b.^2 - 4*a*c))./(2*a); % Solve for omega squared (smallest of the two solutions is of interest here)
% wres = sqrt(wresSquared); % Only positive solution is of interest here
% fres = wres/(2*pi);
% [val ind] = min(abs(fres - f)); % 
% ycdofFull(i) = fres(ind);
% end
% 
% figure;plot(t2,yrb,t2,ycdofFull);legend('RB','CDOF');xlabel('T2');ylabel('wres')
% figure;plot(t2,yrb-ycdofFull);xlabel('T2');ylabel('wresRB-wCDOF')
% 
% 
% 
