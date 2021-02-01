% Script used to load fdComp data and to evaluate power with and without CDOF

%% Clear, clc

clear all; clc;

%% Load WAMIT data

r = 0.88;
mfile_name = sprintf('waveBot_r%.2f_CDOF.mat',r);
load(fullfile('data',mfile_name))

fdComp = fdComp_Rigid;
T = fdComp.T;
f = 1./T;
w = 2*pi*f;

%% Create intrinsic impedance and excitation force matrices/vectors
Zi_rigid = shiftdim(fdComp.B + 1i * ( w .* (fdComp.A + fdComp.M) - (fdComp.C+fdComp.K) ./ w),-2);
Hex_rigid = transpose(fdComp.Fex);

M = shiftdim(repmat(fdComp_CDOF.M,1,1,length(w)),2); C = shiftdim(repmat(fdComp_CDOF.C,1,1,length(w)),2); K = shiftdim(repmat(fdComp_CDOF.K,1,1,length(w)),2);
Zi_cdof = shiftdim(fdComp_CDOF.B + 1i * ( repmat(w,1,fdComp_CDOF.DoF,fdComp_CDOF.DoF) .* (fdComp_CDOF.A + M) - (C+K) ./ repmat(w,1,fdComp_CDOF.DoF,fdComp_CDOF.DoF)),1);
Hex_cdof = transpose(squeeze(fdComp_CDOF.Fex));

% Make test plots for cases with and without the CDOF
figure; grid on; hold on;
pd1 = plot(1./fdComp.T, mag2db(abs(squeeze(Zi_rigid))));
pd2a = plot(1./fdComp_CDOF.T, mag2db(abs(squeeze(Zi_cdof(1,1,:)))),'--');
pd2b = plot(1./fdComp_CDOF.T, mag2db(abs(squeeze(Zi_cdof(1,2,:)))));
pd2c = plot(1./fdComp_CDOF.T, mag2db(abs(squeeze(Zi_cdof(2,2,:)))));
legend([pd1,pd2a,pd2b,pd2c],'Rigid body','CDOF - diag 1','CDOF - off-diag','CDOF - diag 2')
set(gca,'xscale','log')
ylabel('Impedance magnitude [dB]')
xlabel('Frequency [Hz]')

figure; grid on; hold on;
pd1 = plot(1./fdComp.T, abs(fdComp.Fex));
pd2a = plot(1./fdComp.T, abs(squeeze(fdComp_CDOF.Fex(:,1,1))),'--');
pd2b = plot(1./fdComp.T, abs(squeeze(fdComp_CDOF.Fex(:,1,2))));
legend([pd1,pd2a,pd2b],'Rigid body','CDOF - 1','CDOF - 2')
set(gca,'xscale','log')
ylabel('Excitation magnitude [N/m]')
xlabel('Frequency [Hz]')

%% Wave conditions - irregular example
% example_spectra = load('C:\Users\AlfredCotten\Desktop\MATLAB\snl_cyl_array\swan\exampleSpectra.mat');
% 
% S_mean.S = mean(example_spectra.Suni_t,2)/(2*pi);
% S_mean.w = example_spectra.f*2*pi;
% S_mean.type = 'freq';
% 
% figure; grid on; hold on;
% plot(1./f, S_mean.S)
% xlabel('Frequency [Hz]')
% ylabel('Spectral density [?]')
% 
% Te_mean = spec2char(S_mean,5);
% Hm0_mean = spec2char(S_mean,1);
% J_mean = 1025*9.81^2/(64*pi)*Hm0_mean^2 * Te_mean;

for n = 1:length(f) % Optimise power production at each frequency of wave

%% Wave conditions - regular
S_mean.S = zeros(1,length(w)); S_mean.S(n) = 1;
% S_mean.w = 2*pi*0.625; % WaveBot should have natural frequency at 0.625Hz
S_mean.type = 'freq';

H = 2 * S_mean.S(n); % This is just by our chosen convention - this value is used directly as the wave amplitude in the WecPower function.
J(n) = 1025*9.81^2/(32*pi) * T(n) * H^2;

%% Compute power for cases with and without CDOF
[dpow_rigid,~] = runPowStudy(f,Zi_rigid,Hex_rigid,S_mean,'plotFlag',0);
[dpow_cdof,~] = runPowStudy(f,Zi_cdof,Hex_cdof,S_mean,'plotFlag',0);

% fprintf('\nPI Control:\t\t\t\t%f\n')
% i=1;
% fprintf('\n\nDevice power (rigid body):\t\t\t\t%f\n', dpow_rigid(i).P)
% fprintf('Device power (w/ CDOF):\t\t\t\t%f\n', dpow_cdof(i).P)
% fprintf('with CDOF / without CDOF:\t\t\t\t%f\n', dpow_cdof(i).P/dpow_rigid(i).P)
% fprintf('Device power (cc control):\t\t\t\t%f\n', dpow_rigid(i).Pub)
% 
% figure;hold on;
% plot(T,-dpow_rigid(i).P_f)
% plot(T,-dpow_cdof(i).P_f)
% plot(T,-dpow_rigid(i).Pub_f)
% legend('PI control - rigid body','PI control - w CDOF','cc control')
% xlabel('Wave period / s')
% ylabel('Power production / W')
% 
% fprintf('\nP Control:\t\t\t\t%f\n')
% i=2;
% fprintf('\n\nDevice power (rigid body):\t\t\t\t%f\n', dpow_rigid(i).P)
% fprintf('Device power (w/ CDOF):\t\t\t\t%f\n', dpow_cdof(i).P)
% fprintf('with CDOF / without CDOF:\t\t\t\t%f\n', dpow_cdof(i).P/dpow_rigid(i).P)
% fprintf('Device power (cc control):\t\t\t\t%f\n', dpow_rigid(i).Pub)
% 
% plot(T,-dpow_rigid(i).P_f)
% plot(T,-dpow_cdof(i).P_f)
% legend('PI control - rigid body','PI control - w CDOF','cc control','P control - rigid body','P control - w CDOF')

fprintf('\n\nDevice power (cc fraction):\t\t\t\t%f\n', dpow_rigid(2).P./dpow_rigid(2).Pub)

Pow_rigid_PIcntrl(n) = -dpow_rigid(1).P;
Pow_rigid_Pcntrl(n) = -dpow_rigid(2).P;
Pow_rigid_ccCntrl(n) = -dpow_rigid(2).Pub;

Pow_cdof_PIcntrl(n) = -dpow_cdof(1).P;
Pow_cdof_Pcntrl(n) = -dpow_cdof(2).P;
Pow_cdof_ccCntrl(n) = -dpow_cdof(2).Pub;

end

figure;hold on;
plot(1./f,Pow_rigid_PIcntrl,'-r')
plot(1./f,Pow_rigid_Pcntrl,'-b',1./f,Pow_rigid_ccCntrl,'-k');
plot(1./f,Pow_cdof_Pcntrl,'--b',1./f,Pow_cdof_ccCntrl,'--k');
legend('Rigid - PI control','Rigid - P control','Rigid - CC control','CDOF - P control','CDOF - CC control')
title('WaveBot')
xlabel('Wave period / s')
ylabel('Power production / W')
ylim([0 max(Pow_cdof_Pcntrl)*1.5])

% H = 2;
% J = 1025*9.81^2/(32*pi) .* T' * H^2;
% 
figure;hold on;
plot(f,Pow_rigid_Pcntrl./Pow_rigid_ccCntrl,'-b',f,Pow_rigid_ccCntrl./Pow_rigid_ccCntrl,'-k');
plot(f,Pow_cdof_Pcntrl./Pow_rigid_ccCntrl,'--b');
legend('Rigid - P control','Rigid - CC control','CDOF - P control')
title('WaveBot')
xlabel('Wave frequency / Hz')
ylabel('Fraction of Power from CC control')
ylim([0 1.5])

% figure;hold on;
% plot(f,Pow_cdof_Pcntrl,f,Pow_cdof_ccCntrl);
% legend('P control','CC control')
% title('WaveBot with CDOF')
% xlabel('Wave frequency / Hz')
% ylabel('Power production / W')

