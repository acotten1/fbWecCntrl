%{ 
mwave - A water wave and wave energy converter computation package 
Copyright (C) 2014  Cameron McNatt

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Contributors:
    C. McNatt
%}
function [geo] = makePanel_waveBot(r1, r2, t1, t2, Ntheta, Nr, Nz, varargin)

opts = checkOptions({{'Quarter'}, {'NoInt'}}, varargin);
quart = opts(1);
noInt = opts(2);

slope = sqrt((t2-t1)^2+(r2-r1)^2); % Length between two radii, along tapered section

dr(1,1) = r1/Nr(1); dr(2,1) = r2/Nr(2);
r.r1 = r1:-dr(1):0; r.r2 = r2:-dr(2):0;
dzn(1,1) = t1/Nz(1); dzn(2,1) = (t2-t1)/Nz(2);
z.t1 = 0:-dzn(1):-t1; z.t2 = -t1:-dzn(2):-t2;


if (quart)
    if (mod(Ntheta, 4) ~= 0)
        error('The number of panels in the theta direction must be divisible by 4 in order to generate a quarter geometry');
    end
    
    Ntheta = Ntheta/4;
    
    dtheta = pi/2/Ntheta;
    theta = 0:dtheta:(pi/2);
else
    dtheta = 2*pi/Ntheta;
    theta = 0:dtheta:2*pi;
end

% make outer surface
ctheta = cos(theta); ctheta = cos(theta);
stheta = sin(theta); stheta = sin(theta);
cirx.r1 = r1*ctheta; cirx.r2 = r2*ctheta;
ciry.r1 = r1*stheta; ciry.r2 = r2*stheta;

if (z.t1(1) == 0)     % TO DO (may want to remove this section if never panelling up complete wavebot with above water section)
    pans(Ntheta*(sum(Nz) + sum(Nr)),1) = Panel;
    topZero = true;
else    
    if (noInt)
        pans(Ntheta*(sum(Nz) + sum(Nr)),1) = Panel;
    else
        pans(Ntheta*(sum(Nz) + sum(Nr)),1) = Panel; % TO DO?
    end
    topZero = false;
end

np = 0;

for n = 1:Ntheta
    
    % top (interior free surface)
    for m = 1:Nr(2)
        verts = zeros(4,3);
        verts(1,:) = [r.r2(Nr(2)-m+2)*ctheta(n) r.r2(Nr(2)-m+2)*stheta(n) z.t1(1)];
        verts(2,:) = [r.r2(Nr(2)-m+1)*ctheta(n) r.r2(Nr(2)-m+1)*stheta(n) z.t1(1)];
        verts(3,:) = [r.r2(Nr(2)-m+1)*ctheta(n+1) r.r2(Nr(2)-m+1)*stheta(n+1) z.t1(1)];
        verts(4,:) = [r.r2(Nr(2)-m+2)*ctheta(n+1) r.r2(Nr(2)-m+2)*stheta(n+1) z.t1(1)];

        np = np + 1;
        pans(np) = Panel(verts);
        if (topZero)
            pans(np).IsWet = true;
            pans(np).IsInterior = true;
        else
            pans(np).IsWet = false;
        end
        pans(np).IsBody = true;
    end
    
    % bottom
    for m = 1:Nr(1)
        verts = zeros(4,3);
        verts(1,:) = [r.r1(m)*ctheta(n) r.r1(m)*stheta(n) z.t2(end)];
        verts(2,:) = [r.r1(m+1)*ctheta(n) r.r1(m+1)*stheta(n) z.t2(end)];
        verts(3,:) = [r.r1(m+1)*ctheta(n+1) r.r1(m+1)*stheta(n+1) z.t2(end)];
        verts(4,:) = [r.r1(m)*ctheta(n+1) r.r1(m)*stheta(n+1) z.t2(end)];

        np = np + 1;
        pans(np) = Panel(verts);
        pans(np).IsWet = true;
        pans(np).IsInterior = false;
        pans(np).IsBody = true;
    end
    
    % top wall
    for m = 1:Nz(1)
        verts = zeros(4,3);
        verts(1,:) = [cirx.r2(n) ciry.r2(n) z.t1(m)];
        verts(2,:) = [cirx.r2(n) ciry.r2(n) z.t1(m+1)];
        verts(3,:) = [cirx.r2(n+1) ciry.r2(n+1) z.t1(m+1)];
        verts(4,:) = [cirx.r2(n+1) ciry.r2(n+1) z.t1(m)];

        np = np + 1;
        pans(np) = Panel(verts);
        if (z.t1(m) <= 0)
            pans(np).IsWet = true;
        else
            pans(np).IsWet = false;
        end
        pans(np).IsInterior = false;
        pans(np).IsBody = true;
    end
    
    % lower, sloped/tapered wall
    for m = 1:Nz(2)
        r_temp = r2 - ([m-1;m]/Nz(2)*(r2-r1));
        cirx_temp = r_temp*ctheta; % the radius is now dependent on height for this surface
        ciry_temp = r_temp*stheta;
        verts = zeros(4,3);
        verts(1,:) = [cirx_temp(1,n) ciry_temp(1,n) z.t2(m)];
        verts(2,:) = [cirx_temp(2,n) ciry_temp(2,n) z.t2(m+1)];
        verts(3,:) = [cirx_temp(2,n+1) ciry_temp(2,n+1) z.t2(m+1)];
        verts(4,:) = [cirx_temp(1,n+1) ciry_temp(1,n+1) z.t2(m)];

        np = np + 1;
        pans(np) = Panel(verts);
        if (z.t2(m) <= 0)
            pans(np).IsWet = true;
        else
            pans(np).IsWet = false;
        end
        pans(np).IsInterior = false;
        pans(np).IsBody = true;
    end
    
end

if (quart)
    geo = PanelGeo(pans, 'Xsym', 'Ysym');
else
    geo = PanelGeo(pans);
end


end