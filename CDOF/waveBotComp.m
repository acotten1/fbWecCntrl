function [fdComp, hydroBody, waveField] = waveBotComp(rho, h, T, r1, r2, t1, t2, t3, varargin)
% WAVEBOTCOMP [fdComp, hydroBody, waveField] = waveBotComp(rho, h, T, 
%       r1, r2, t1, t2, t3, varargin)
%
%   Inputs:
%       rho: water density
%       h: water depth (cannot be Inf)
%       T: vector of wave periods
%       The following four parameters describe the submerged portion 
%       of the WaveBot body:
%       r1 - smaller (base) radius
%       r2 - larger (top) radius
%       t1 - Height of wide cylindrical part
%       t2 - Total height of submerged part (t1 + height of tapered part)
%   
%       t3 is the height of the cylindrical part above water
%       varargin: optional arguments:
%           'motions', (1,6): 6Dof motions vector: e.g. [1 1 1 1 1 1] -
%               motions in all 6 Dof. [0 0 1 0 0 0] - motion in only heave
%           'panelSize', (1,1): target panel size of the WAMIT geometry
%           'hydroBody': compute the hydroBody
%           'waveField', (1,2), (1,2): compute the wave field. Inpute args
%               are the [length in x, length in y] and [N pnt x, N pnt y]
%           'array', (N,2): compute an array of cylinders (directly with
%               WAMIT). Input arg is XY locations of cylinders in array
%           'wamitPath', (1,:): user defined WAMIT .exe path. Input arg is
%               path string. Default: 'C:\wamitv7'
%           'directSolver': use the Wamit direct solver, rather than the
%               iterative solver. The direct solver takes longer than the 
%               iterative, but is guanteed to find a solution. 
%
%   Outputs:
%       fdComp: frequency domain computation. FreqDomComp object
%       hydroBody: "hydroBody". HydroBody object. Used to compute arrays
%           with FreqDomArrayComp. If optional arg 'hydroBody' is not 
%           present, then is empty
%       waveField: the wave field computed by WAMIT. FBWaveField object. If
%           optional arg 'waveField' is not present, then is empty
%
%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%
%   Computes for the WaveBot geometry as defined by the dimensions 
%   specified above:
%       - linear wave hydrodynamic coefficients using WAMIT
%       - mass properties assuming constant density, weight = displacement
%       - so-called "hydroBody" which is an object that contains the
%       "diffraction transfer matrix" and "radiation coefficients", that
%       can be used to compute array interactions using interaction theory
%       - wave fields using WAMIT
%       - arrays of the cylinder using the so-called direct method. i.e.
%       using WAMIT to compute linear wave coefficients for all bodies in
%       array.
%
%   Notes:
%       - the default size of the panels and the spacing of the points used
%       to compute the hydroBody is driven by the shortest wavelength. 
%       - this function will create a subfolder called 'wamrun' in the
%       directory where this files is located and put the WAMIT i/o files
%       there.
%
%   References:
%       McNatt, J. C., V. Venugopal, D. Forehand. "A novel method for 
%           deriving the diffraction transfer matrix and its application to
%           multi-body interactions in water waves". Ocean Engineering.
%           2015. 94, 173-185
%       McNatt, J. C. "Cylindrical Linear Water Waves and their Application 
%           to the Wave-body Problem". 2016. PhD thesis. University of
%           Edinburgh. https://era.ed.ac.uk/handle/1842/20378
%
%   See Also: \mwave\Examples
%
%   Author: C. McNatt. cameron.mcnatt@mocean.energy
%   Updated: June 2020

% Check input arguments
arguments
    rho(1,1) double {mustBePositive, mustBeFinite}
    h(1,1) double {mustBePositive, mustBeFinite}
    T(1,:) double {mustBePositive, mustBeFinite}
    r1(1,1) double {mustBePositive, mustBeFinite}
    r2(1,1) double {mustBePositive, mustBeFinite}
    t1(1,1) double {mustBePositive, mustBeFinite}
    t2(1,1) double {mustBePositive, mustBeFinite}
    t3(1,1) double {mustBePositive, mustBeFinite}
end
arguments (Repeating)
    varargin
end

[opts, args] = checkOptions({{'motions', 1}, {'panelSize', 1}, ...
    {'hydroBody'}, {'waveField', 2}, {'array', 1}, {'wamitPath', 1}, ...
    {'directSolver'},{'cdofVol',1}}, varargin);

mots = [1 1 1 1 1 1];
if opts(1)
    mots = args{1};
    validateattributes(mots, {'numeric'}, {'binary'});
    if length(mots)~=6 && length(mots)~=7
        error('motions vector should be either 6 or 7 in length.')
    end
end
panSize = [];
if opts(2)
    panSize = args{2};
    validateattributes(panSize, {'numeric'}, {'size', [1 1], '>', 0})
end
compHb = opts(3);
compWf = opts(4);
if compWf
    lenWf = args{4}{1};
    Nwf = args{4}{2};
    validateattributes(lenWf, {'numeric'}, {'size', [1 2], '>', 0})
    validateattributes(Nwf, {'numeric'}, {'size', [1 2], '>', 0, 'integer'})
end
compArray = opts(5);
if compArray 
    arrayXy = args{5};
    validateattributes(arrayXy, {'numeric'}, {'ncols', 2})
end
wamPath = [];
if opts(6)
    wamPath = args{6};
    validateattributes(wamPath, {'char'}, {'nonempty'})
end
directSolver = opts(7);

if opts(8)
    cdofVol = args{8};
else
    cdofVol = 0.01; % Set close to zero to approximate the rigid body case by default.
end

if compHb && compArray
    error('cannot compute ''hydroBody'' and ''array'' options at the same time.');
end  

% Set fundamental constants for use in computing the extra stiffness of the compressible DOF 
pAtm = 101325; % Atmospheric pressure(in Pascals)
adInd = 1.4; % Adiabatic index for air
g = 9.80665; % gravitational acceleration

% Set up cylinder geometry

% Get the shortest wavelength to determine default value panel size and
% wave field point spacing values
lam = min(IWaves.T2Lam(T, h));
delta = min([lam/4, r1/4, t2/5, 2*pi*r2/8]);
% delta = 0.03;
if isempty(panSize)
    panSize = delta;
end

Ntheta = round(2*pi*r1/panSize); % For simplicity of meshing the sloped surface, take Ntheta from the smaller radius.
Nr(1,1) = round(r1/panSize);
Nr(2,1) = round(r2/panSize);
Nz(1,1) = round(t1/panSize);
Nz(2,1) = round(sqrt((t2-t1)^2+(r2-r1)^2)/panSize);

wec = FloatingWaveBot(rho, r1, r2, t1, t2, t3, Ntheta, Nr, Nz); % Add 'UseSym' option here to only panel a quarter of the geometry in the WAMIT .GDF files.
wec.Modes = ModesOfMotion(mots);   % compute for heave only
wec.Dpar = zeros(size(wec.Dpto)); % Dpto, K and M are reset for total no. of modes in the above line 
                        % by calling the set.Modes method within the FloatingBody class, but Dpar isn't.
if wec.Ngen > 0 % i.e. if using CDOF
    wec.K(7,7) = (pi*r1^2)^2*adInd*(pAtm + rho*g*t2)/cdofVol; % stiffness term computed for the compressible volume. See Kurniawan paper or SAND report for derivation. 
end

% Set up wamit computation
runName = 'waveBotComp';
folder = [fileparts(which('waveBotComp')) '\wamrun'];
if 0 == exist(folder, 'dir')
    mkdir(folder);
end
delete([folder '\*']);

wamRun = WamitRunCondition(folder, runName);  

if ~isempty(wamPath)
    wamRun.ExePath = wamPath;
    wamRun.ScratchPath = [wamPath '\scratch']; 
    wamRun.UseridPath = wamPath;        % Location of UserId (license)
end
wamRun.UseDirectSolver = directSolver;

wamRun.Rho = rho;      
wamRun.T = T;                
wamRun.Beta = 0;              
wamRun.H = h;    
wamRun.NCPU = 6;
wamRun.RAMGBmax = 8;

if compArray
    wec0 = FloatingBody(wec);
    Nb = size(arrayXy, 1);
    
    wec(Nb, 1) = FloatingBody;
    for n = 1:Nb
        wec(n) = FloatingBody(wec0);
        wec(n).XYpos = arrayXy(n, :);
    end
end
wamRun.FloatingBodies = wec; 

if compHb
    r = wec.Rcir;                   
    dr = delta;                       
    r = r + dr;                     

    nZ = round(h/(delta/2));                 
    nTheta = round(2*pi*r/delta);   % was nTheta = 2^8;  
    nTheta = 2^ceil(log(nTheta)/log(2));

    cylArray = BemCylArray(r, nTheta, nZ);            
    wamRun.CylArray = cylArray;   
end

if compWf
    start = [-lenWf(1)/2 -lenWf(2)/2 0];    % start location, [x, y, z]
    delta = [lenWf./(Nwf - 1) 1];           % delta (point spacing) in [x, y, z]
    nPoint = [Nwf(1) Nwf(2) 1];             % Number of points in [x, y, z]. Having only 1  
                                            % point in the z direction creates a 2D array
    wamRun.FieldArray = BemFieldArray(start, delta, nPoint);
end

if compArray
end
   
wamRun.WriteRun;

if wec.Ngen > 0 && wec.WamIGenMds == 0 % i.e. if using defmod approach with some generalised modes
    writeDefmod(t2); % Write and compile defmod based on the current geometry 
              % (input argument should be the z-location of the panels associated with the heave DOF).
    overwriteBatch(); % Modify the batch file so gdf filename is automatically provided
end

% Run Wamit
wamRun.Run;                                      

% Read results, and save useful objects
wamResult = WamitResult(wamRun);  
wamResult.ReadResult;          


forces = wamResult.FreqDomForces;
fdComp = FreqDomComp(forces, wec);

hydroBody = [];
if compHb
    waveCir = wamResult.WavePoints;
    % create HydroBody and save it
    hydroBody = computeHydroBody(waveCir, forces, wec, 'AxisSym', wec.Modes, ...
        'SigFigCutoff', 5, 'AccTrim');
end

waveField = [];
if compWf
    waveField = wamResult.WaveArray;  
end
