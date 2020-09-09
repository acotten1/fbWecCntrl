% Tabulates performance for complex-conjugate (CC - perfect impedance
% matching) and proportional-integral (PI) controllers designed based on
% either the hydro-mechanical system (e.g., "CC on mech") or the electrical
% system (e.g., "CC on elec") in a single sea state.

% clc
clear
close all

optimOpts = optimoptions('fminunc',...
    'MaxFunctionEvaluations',1e6, 'MaxIterations', 1e6, 'Display', 'off');

%% Load WEC device data

cf = 60;
mf = load('waveBot_heaveModel.mat');
Zi = mf.Zi_frf(cf:end,1);
Hex = mf.H_frf(cf:end,1)*1e1;
f = mf.f(cf:end,1);
w = 2*pi*f;
dw = w(2)-w(1);

Zpto = PTO_Impedance(w,[1, 0, 0, 0, sqrt(2/3), 1e-3, 0]); % [N, Id, Bd, Kd, Kt, Rw, Lw]

%% Define sea state and excitation

Hs = 0.125;
Tp = 2;
gamma = 3.3;

S = jonswap(w, [Hs, Tp, gamma]);    % Wave energy density spectrum
A = sqrt(2*dw*S.S(:));              % wave amplitude spectrum
Fe = A .* Hex(:);

Pmax = abs(Fe).^2 ./ (8*real(Zi));

%% Design controllers

%---------------------------------
wc(1).leg = 'CC on mech';
wc(1).ZL = Zi2ZL(Zpto, conj(Zi));

%---------------------------------
wc(2).leg = 'PI on mech';
wc(2).cinfo.type = 'PI';
wc(2).cinfo.w = w;
wc(2).cinfo.x0 = ones(1,2)*0.1;
wc(2).objfun = @(x) Pmech( Zi2ZL(Zpto,fbc(x,wc(2).cinfo)),...
    Zpto,...
    Zi,Fe );
[wc(2).y, wc(2).fval] = fminunc(wc(2).objfun, wc(2).cinfo.x0, optimOpts);
wc(2).ZL = Zi2ZL(Zpto,fbc(wc(2).y, wc(2).cinfo));

%---------------------------------
wc(3).leg = 'CC on elec';
wc(3).ZL = conj( squeeze(Zpto(2,2,:)) ...
    - squeeze(Zpto(1,2,:)) .* squeeze(Zpto(2,1,:)) ...
    ./ (squeeze(Zpto(1,1,:)) + Zi) );

%---------------------------------
wc(4).leg = 'PI on elec';
wc(4).cinfo.type = 'PI';

wc(4).cinfo.w = w;
wc(4).cinfo.x0 = ones(1,2);
wc(4).objfun = @(x) Pelec( Zi2ZL(Zpto,fbc(x,wc(4).cinfo)),...
    Zpto,...
    Zi,Fe );
[wc(4).y, wc(4).fval] = fminunc(wc(4).objfun, wc(4).cinfo.x0, optimOpts);
wc(4).ZL = Zi2ZL(Zpto,fbc(wc(4).y, wc(4).cinfo));


%% Calculate power and efficiency

for ii = 1:length(wc)
    
    % evaluate mech performance
    [Pmech_tot(ii), mPmech(:,ii)] = Pmech(wc(ii).ZL, Zpto, Zi, Fe);
    assert(-1*Pmech_tot(ii) <= sum(Pmax),...
        sprintf('''%s'' making more mechanical power than theoretical limit',wc(1).leg))
    
    % evaluate mech performance
    [Pelec_tot(ii), mPelec(:,ii)] = Pelec(wc(ii).ZL, Zpto, Zi, Fe);
    assert(-1*Pelec_tot(ii) <= sum(Pmax),...
        sprintf('''%s'' making more mechanical power than theoretical limit',wc(1).leg))
    
    legCel{ii} = wc(ii).leg;
end

eta_mech = Pmech_tot./(-1 * sum(Pmax));
eta_elec = Pelec_tot./(-1 * sum(Pmax));

%% Tabulate results

T = table(-1*Pmech_tot(:)/1e3,eta_mech(:),-1*Pelec_tot(:)/1e3,eta_elec(:),...
    'VariableNames',{'MechPow_kW','MechEfficiency','ElecPow_kW','ElecEffciency'},...
    'RowNames',legCel);
disp(T)

clear input
input.data = T{:,1:end};
input.dataFormat = {'%.2g',};
input.tableColumnAlignment = 'c';
input.tableBorders = 0;
input.makeCompleteLatexDocument = 0;
input.tableRowLabels = T.Properties.RowNames;
input.tableCaption = 'Text caption';
input.tableLabel = 'tab:TestLabel';
input.tableColLabels = T.Properties.VariableNames;
latex = latexTable(input);

%% Plot results

figure('name','Frequency dependent power')
hold on
grid on
plot(mPmech,'-')
ax = gca;
ax.ColorOrderIndex = 1;
plot(mPelec,'--')

figure('name','Bode plot')
hold on
grid on
opt = bodeoptions;
opt.FreqUnits = 'Hz';
opt.Grid = 'on';
for ii = 1:length(wc)
    bodeplot(frd(wc(ii).ZL,w),opt)
    legCel{ii} = wc(ii).leg;
end
legend(legCel)
