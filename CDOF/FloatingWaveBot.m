classdef FloatingWaveBot < FloatingBody
    
    properties (Access = protected)
        r1;
        r2;
        t1;
        t2;
        t3;
    end

    properties (Dependent)
        R1;
        R2;
        T1;
        T2;
        T3;
    end
    
    methods
        % Constructor 
        function [fb] = FloatingWaveBot(rho, r1, r2, t1, t2, t3, varargin)
            % Floating WaveBot of uniform density that rotates about its
            % cg
            % The following four parameters describe the submerged portion 
            % of the WaveBot body.
            % r1 - smaller (base) radius
            % r2 - larger (top) radius
            % t1 - Height of wide cylindrical part
            % t2 - Total height (t1 + height of tapered part)
            % 
            % t3 is the height of the cylindrical part above water
            % optional arguments - number or panels in theta (Ntheta), r
            % (Nr), and z (Nz)
            
            fb = fb@FloatingBody();
            
            % Compute Cg assuming mass is evenly
            % distributed through total WaveBot volume - i.e. compute
            % centre of volume
            vol_cyl = pi*(r2^2*(t1+t3)); % Volume of total cylindrical part (both above and below water)
            vol_truncCone = pi*((1/3)*(r2^2 + r1*r2 + r1^2)*(t2-t1));
            vol_total = vol_cyl + vol_truncCone;
            zcg0_cyl = -t1 + (t1+t3)/2;
            zcg0_truncCone = -t1 -(t2-t1)*r1/(r1+r2);
            zcg0 = zcg0_cyl*(vol_cyl/vol_total) + zcg0_truncCone*(vol_truncCone/vol_total); % Total cg in z-direction, relative to waterline
            
            fb.cg = [0 0 zcg0];

            fb.m = FloatingWaveBot.MassMatrix(rho, r1, r2, t1, t2); % Note: inertia terms are currently set = 0: only use translational DoFs.
            fb.position = [0 0 0];
            
            fb.r1 = r1;
            fb.r2 = r2;
            fb.t1 = t1;
            fb.t2 = t2;
            
            if (~isempty(varargin))
                if (length(varargin) < 3)
                    error('There must be three optional inputs to FloatingCylinder: Ntheta, Nr, Nz');
                end
                Ntheta = varargin{1};
                Nr = varargin{2};
                Nz = varargin{3};
                
                opts = checkOptions({{'UseSym'}, {'NoInt'}}, varargin);
                optsIn = {};
                n = 0;
                if (opts(1))
                    n = n + 1;
                    optsIn{n} = 'Quarter';
                end
                if (opts(2))
                    n = n + 1;
                    optsIn{n} = 'NoInt';
                    fb.iSurfPan = 0;
                else
                    fb.iSurfPan = 1;
                end
                
                panGeo = makePanel_waveBot(fb.r1, fb.r2, fb.t1, fb.t2, Ntheta, Nr, Nz, optsIn{:});
                fb.panelGeo = panGeo;
                fb.iLowHi = 0;
            end
            
            pts = CirContSurf.Compute(fb.cg, r2, 30);
            pts = [pts; r2 0 0];
            fb.wpSec = pts(:, 1:2);
        end
    end
    
    methods
        function [r] = get.R1(fb)
            r = fb.r1;
        end
        
        function [r] = get.R2(fb)
            r = fb.r2;
        end
        
        function [h] = get.T1(fb)
            h = fb.t1;
        end
        
        function [h] = get.T2(fb)
            h = fb.t2;
        end
        
        function [h] = get.T3(fb)
            h = fb.t3;
        end
        
        function [] = MakePanelGeometry(fb, Ntheta, Nr, Nz, varargin)
            opts = checkOptions({'UseSym'}, varargin);
            
            if (opts(1))
                panGeo = makePanel_waveBot(fb.r1, fb.r2, fb.t1, fb.t2, Ntheta, Nr, Nz, 'Quarter');
            else
                panGeo = makePanel_waveBot(fb.r1, fb.r2, fb.t1, fb.t2, Ntheta, Nr, Nz);
            end
            %panGeo.Translate(-fb.position);
            fb.panelGeo = panGeo;
            fb.iLowHi = 0;
        end
    end
    
    methods (Access = protected)
%         function [] = onModifyCg(fb, cg)
%             mass = fb.m(1);
%             wetVol = pi*fb.radius^2*fb.draft;
%             rho = mass/wetVol;
%             zcg = cg(3);
% 
%             fb.m = FloatingCylinder.MassMatrix(rho, fb.radius, fb.height, fb.draft, zcg);
%         end
    end
    
    methods (Static)
        function [M] = MassMatrix(rho, r1, r2, t1, t2, varargin)
            if (~isempty(varargin))
                % zcg is the vertical position in body coordinates, 
                % where z=0 is the calm water free surface
                zcg0 = -draft+height/2;
                zcg = varargin{1};
                delzcg = zcg - zcg0;
            else
                delzcg = 0;
            end
            % rho is the fluid density
            % body mass equals buoyant force
            wetVol = pi*(r2^2*t1 + (1/3)*(r2^2 + r1*r2 + r1^2)*(t2-t1)); % Vol of cylinder + vol of truncated cone
            m = rho*wetVol; % Total mass (Archimedes' principle)

            Izz = 0; % For now, only set up to look at heave motions, so okay to have these set to zero
            Ixx = 0;
            Iyy = Ixx;
            
            if (delzcg ~= 0)
                Ixx = Ixx + m*delzcg^2;
                Iyy = Iyy + m*delzcg^2;
            end

            M = zeros(6,6);

            M(1,1) = m;
            M(2,2) = m;
            M(3,3) = m;
            M(4,4) = Ixx;
            M(5,5) = Iyy;
            M(6,6) = Izz;
        end
    end
end